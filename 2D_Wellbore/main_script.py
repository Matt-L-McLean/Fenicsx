#-------------------------------------
# This script solves the plane-strain problem
# of a wellbore in a poroelastic medium
# Units: stress-Pa, length-m, time-sec.
# Sign convention: compressive stress & contraction strain are positive
#-------------------------------------

# Mesh from Gmsh
#-------------------------------------
from mesh_creation import Gmsh_model
ms,ms_w,xl,yl,dw = 3.5, 0.015, 25, 25, 0.2032 # mesh size at corners [m], well mesh size [m], x-half-length [m], y-half-length [m], wellbore diameter [m]
model,gdim = Gmsh_model(ms,ms_w,xl,yl,dw)
from mesh_conversion import msh_to_xdmf
mesh,cell_tags = msh_to_xdmf(model, gdim) 
from inline_plot import plot_mesh
plot_mesh(mesh, cell_tags)

# Import several Fenicsx/dolfinx packages
#-------------------------------------
from dolfinx import plot, io
from dolfinx.fem import (dirichletbc, Expression, Function, FunctionSpace, 
                         VectorFunctionSpace, TensorFunctionSpace, locate_dofs_topological, Constant)
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.mesh import locate_entities_boundary, locate_entities, meshtags
from ufl import (TestFunction, TrialFunction, dot, dx, grad, inner, nabla_div, div, Identity, sym, Measure,
                 SpatialCoordinate, lhs, rhs, as_vector, tr, FiniteElement, VectorElement, MixedElement, 
                 split, FacetNormal, TensorElement, as_matrix, as_tensor, atan_2, cos, sin)
from petsc4py.PETSc import ScalarType, Options
import numpy as np
import dolfinx
from mpi4py import MPI

# Time stepping properties (Euler method)
#-------------------------------------
t = 0.0 # initialize time [s]
tmax = 86400 # 1 day = 86,400 seconds (units of time are seconds!)
num_steps = 25 # number of steps
dt = tmax/num_steps # time increment

# Defined poroelastic constants
#-------------------------------------
E = 10e9 # Young's modulus [Pa]
nu = 0.22 # Poisson's ratio [-]
phi = 0.2 # Porosity [-]
k_ = 1e-16 # permeability [m^2]
mu_ = 1e-3 # Viscosity [Pa*s]
Ks = 37e9 # Solid bulk modulus [Pa]
Kf = 2.15e9 # Fluid bulk modulus [Pa]

# Derived poroelastic constants
#-------------------------------------
lambda_ = E*nu/((1+nu)*(1-2*nu)) # Lame parameter [Pa]
mu = E/(2*(1+nu)) # Shear modulus [Pa]
model = "plane_strain"
if model == "plane_stress":
    lambda_ = 2*mu*lambda_/(lambda_+2*mu) # Lame parameter for plane-stress elasticity [Pa]
kappa = k_/mu_ # Bulk fluid mobility [m^2/Pa/s]
Km = E/(3*(1-2*nu)) # Drained bulk modulus [Pa]
alpha = 1-Km/Ks # Biot coefficient [-]
M = (phi/Kf + (alpha-phi)/Ks)**(-1) # Biot modulus [Pa]
Ku = Km + alpha**2 * M # Undrained bulk modulus [Pa]
nu_u = (3*Ku-2*mu)/(2*(3*Ku+mu)) # Undrained Poisson's ratio [-]

# FEM space
#-------------------------------------
disp = VectorElement("CG", mesh.ufl_cell(), 2) # Vector piecewise quadratic Lagrange element for u=(u,v)
pres = FiniteElement("CG", mesh.ufl_cell(), 1) # Scalar piecewise linear Lagrange element for p
ten = TensorElement("CG", mesh.ufl_cell(), 2) # Tensor piecewise quadratic Lagrange element for S=(Sxx, Sxy, Sxz,
#                                                                                                  Syx, Syy, Syx,
#                                                                                                  Szx, Szy, Szz)
V0 = MixedElement([disp, pres]) # Solution vector is ordered sol=(u,v,p)
V0 = FunctionSpace(mesh,V0) # Mixed function space
S = FunctionSpace(mesh, ten) # Tensor function space
C = FunctionSpace(mesh, pres) # Scalar function space

# In-situ stress (stress IC *NOT* displacement IC)
# This is added to the stress tensor
# See Coussy,2004 Poromechanics text:
# sigma - sigma_i = C:eps - alpha*(p-p_i)
# Hence, this is NOT a body force but rater an IC
#-------------------------------------
Sv = 35e6 # Total vertical stress [Pa]
SHmax = 29e6 # Total maximum horizontal stress [Pa]
Shmin = 21e6 # Total minimum horizontal stress [Pa]
Pp = 15e6 # Pore pressure [Pa]
Pw = 20e6 # Wellbore pressure [Pa]

Sxx_i = Constant(mesh,Shmin) # x-dir ---> direction of Shmin
Sxy_i = Constant(mesh,0.0) # No initial shear stress
Syy_i = Constant(mesh,SHmax) # y-dir ---> direction of Shmax
p_i = Constant(mesh,Pp) # Initial pore pressure
IC = as_tensor([[Sxx_i, Sxy_i],[Sxy_i, Syy_i]]) # 2D initial stress tensor for S=(Sxx, Sxy,
#                                                                                 Syx, Syy)

# Stress and strain definition
#-------------------------------------
def epsilon(u):
    e = sym(grad(u)) # Strain = 0.5*(grad(u) + grad(u).T)
    return as_tensor([[e[0, 0], e[0, 1]],
                      [e[0, 1], e[1, 1]]])
def sigma(u,p):
    sig = -lambda_*nabla_div(u)*Identity(u.geometric_dimension()) - 2*mu*epsilon(u) + alpha*(p-p_i)*Identity(u.geometric_dimension()) # Total stress tensor
    return as_tensor([[sig[0, 0], sig[0, 1]],
                      [sig[0, 1], sig[1, 1]]]) + IC
def sigma_eff(u,p):
    sig_e = -lambda_*nabla_div(u)*Identity(u.geometric_dimension()) - 2*mu*epsilon(u) + alpha*(p-p_i)*Identity(u.geometric_dimension()) - p*Identity(u.geometric_dimension()) # Effective stress tensor
    return as_tensor([[sig_e[0, 0], sig_e[0, 1]],
                      [sig_e[0, 1], sig_e[1, 1]]]) + IC
def epsilon_v(u):
    return nabla_div(u) # Volumetric strain

# Split mixed FEM space
#-------------------------------------
V = TestFunction(V0) # Test function of mixed element space
del_u, del_p = split(V) # Test variables for disp. and pressure

# We now need a way to keep track of current and previous solutions
# to discretize time derivatives in fluid mass balance

wk1 = Function(V0) # Current solution vector
wk = Function(V0) # Previous solutiun vector

uk1, pk1 = split(wk1) # Split into seperate variables
uk, pk, = split(wk) # Split into seperate variables

# Initial Conditions --> displacement and pressure *ONLY*
# Remember, stress is already initialized
# But, initial pore pressure must be defined
#-------------------------------------
wk1.x.array[:] = 0.0 # Initialize 1st time step
wk1.sub(0).sub(0).interpolate(lambda x: np.zeros(x.shape[1])) # Access to x-displacement IC
wk1.sub(0).sub(1).interpolate(lambda x: np.zeros(x.shape[1])) # Access to y-displacement IC
wk1.sub(1).interpolate(lambda x: Pp*np.ones(x.shape[1])) # Access to pressure IC
wk1.x.scatter_forward()

# Solve linear variational problem
#-------------------------------------
F = kappa*dot(grad(pk1),grad(del_p))*dx + alpha*((nabla_div(uk1)-nabla_div(uk))/dt)*del_p*dx + alpha**2/(Ku-Km)*(pk1-pk)/dt*del_p*dx + inner(sigma(uk1,pk1),epsilon(del_u))*dx # Residual

# Boundary conditions
#-------------------------------------
x = SpatialCoordinate(mesh) # Locate x,y coordinates of elements
n = -FacetNormal(mesh) # normal direction to elements
norm = as_vector([n[0], n[1]]) # This is needed for Neumann BCs
theta = atan_2(x[1],x[0]) # Radial coordinate "theta"
rotate_1 = as_tensor([[cos(theta), sin(theta)], # Rotation matrix cartesion --> radial
                      [-sin(theta), cos(theta)]])
rotate_2 = as_tensor([[cos(theta), -sin(theta)], # Rotation matrix cartesion --> radial
                      [sin(theta), cos(theta)]])

boundaries = [(1,lambda x: np.isclose(x[0], -xl)),
                  (2,lambda x: np.isclose(x[1], -yl)),
                  (3,lambda x: np.isclose(x[0], xl)),
                  (4,lambda x: np.isclose(x[1], yl)),
                  (5,lambda x: np.isclose(np.sqrt(x[0]**2 + x[1]**2),dw/2))] # pick out the left, bottom, right, top, and wellbore boundary DOFs

facet_indices, facet_markers = [], []
fdim = mesh.topology.dim - 1
for (marker, locator) in boundaries:
    facets = locate_entities(mesh, fdim, locator)
    facet_indices.append(facets)
    facet_markers.append(np.full(len(facets), marker))
facet_indices = np.array(np.hstack(facet_indices), dtype=np.int32)
facet_markers = np.array(np.hstack(facet_markers), dtype=np.int32)
sorted_facets = np.argsort(facet_indices)
facet_tag = meshtags(mesh, fdim, facet_indices[sorted_facets], facet_markers[sorted_facets])

ds = Measure("ds", domain=mesh, subdomain_data=facet_tag) # Integration measure

class BoundaryCondition():
    def __init__(self, type, marker, values, ind_1, ind_2):
        self._type = type
        if type == "Dirichlet":
            u_D = values
            facets = np.array(facet_tag.indices[facet_tag.values == marker])
            if ind_1 == 0:
                dofs = locate_dofs_topological(V0.sub(ind_1).sub(ind_2), fdim, facets)
                self._bc = dirichletbc(u_D, dofs, V0.sub(ind_1).sub(ind_2))
            else:
                dofs = locate_dofs_topological(V0.sub(ind_1), fdim, facets)
                self._bc = dirichletbc(u_D, dofs, V0.sub(ind_1))
        elif type == "Neumann":
            if ind_1 == 0:
                self._bc = values*inner(norm, del_u) * ds(marker)
            else:
                self._bc = inner(values, del_p) * ds(marker)
        else:
            raise TypeError("Unknown boundary condition: {0:s}".format(type))
    @property
    def bc(self):
        return self._bc

    @property
    def type(self):
        return self._type

# Fixed displacement along outer boundaries
# Mechanical load (from drilling fluid) along wellbore
# Constant pore pressure along wellbore
# Constant (initial) pore pressure along outer boundaries
# Could change last BC to zero flux for shales

boundary_conditions = [BoundaryCondition("Dirichlet", 1, ScalarType(0), 0, 0),
                           BoundaryCondition("Dirichlet", 2, ScalarType(0), 0, 1),
                           BoundaryCondition("Dirichlet", 3, ScalarType(0), 0, 0),
                           BoundaryCondition("Dirichlet", 4, ScalarType(0), 0, 1),
                           BoundaryCondition("Neumann", 5, ScalarType(Pw), 0,0),
                           BoundaryCondition("Dirichlet", 5, ScalarType(Pw), 1, 0),
                           BoundaryCondition("Dirichlet", 1, ScalarType(Pp), 1, 0),
                           BoundaryCondition("Dirichlet", 2, ScalarType(Pp), 1, 0),
                           BoundaryCondition("Dirichlet", 3, ScalarType(Pp), 1, 0),
                           BoundaryCondition("Dirichlet", 4, ScalarType(Pp), 1, 0)]

bcs = []
for condition in boundary_conditions:
    if condition.type == "Dirichlet":
        bcs.append(condition.bc)
    else:
        F += condition.bc


# Transient solution
#-------------------------------------
file = io.XDMFFile(mesh.comm, "Solution.xdmf", "w") # Output file --> visualize in Paraview
file.write_mesh(mesh) # Attach mesh to the file

uk1_sol = wk1.sub(0) # Current displacement solution vector --> t = k+1
pk1_sol = wk1.sub(1) # Current pressure solution vector --> t = k+1

uk1_sol.name = "Displacement"
pk1_sol.name = "Pore Pressure"

wk.x.array[:] = wk1.x.array # Initialize previous solution

problem = NonlinearProblem(F, wk1, bcs) # Newton iteration in Dolfinx requires nonlinear problem
solver = NewtonSolver(mesh.comm, problem) # Newton method
solver.convergence_criterion = "incremental" # Use either "incremental" or "residual"
solver.rtol = 1e-11 # Relative tolerance
solver.report = True # Output solver report
solver.max_it = 50 # Stop if no solution within 50 iterations

ksp = solver.krylov_solver
opts = Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "preonly"
opts[f"{option_prefix}pc_type"] = "lu"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()
        
while (t<tmax):
    t += dt # Update time
    r = solver.solve(wk1) # Solve the problem!
    
    print(f"Step {int(t/dt)}: num iterations: {r[0]}")
    
    u_sol,p_sol = split(wk1) # Nodal solution at t = k+1
    
    # Output stress and stain --> visualize in Paraview
    expr_1 = Expression(sigma(u_sol,p_sol), S.element.interpolation_points)
    stress = Function(S)
    stress.interpolate(expr_1)
    stress.name = "Total Stress"
    
    expr_2 = Expression(sigma_eff(u_sol,p_sol), S.element.interpolation_points)
    stress_eff = Function(S)
    stress_eff.interpolate(expr_2)
    stress_eff.name = "Effective Stress"
    
    expr_3 = Expression(epsilon(u_sol), S.element.interpolation_points)
    strain = Function(S)
    strain.interpolate(expr_3)
    strain.name = "Strain"
    
    expr_4 = Expression(epsilon_v(u_sol), C.element.interpolation_points)
    strain_v = Function(C)
    strain_v.interpolate(expr_4)
    strain_v.name = "Volumetric Strain"
    
    expr_5 = Expression(rotate_1*sigma(u_sol,p_sol)*rotate_2, S.element.interpolation_points)
    stress_r = Function(S)
    stress_r.interpolate(expr_5)
    stress_r.name = "Polar Total Stress"
    
    expr_6 = Expression(rotate_1*epsilon(u_sol)*rotate_2, S.element.interpolation_points)
    strain_r = Function(S)
    strain_r.interpolate(expr_6)
    strain_r.name = "Polar Strain"
    
    # Write all of the above to output file
    file.write_function(uk1_sol,t)
    file.write_function(pk1_sol,t)
    file.write_function(stress,t)
    file.write_function(stress_eff,t)
    file.write_function(strain,t)
    file.write_function(strain_v,t)
    file.write_function(stress_r,t)
    file.write_function(strain_r,t)
    
    wk.x.array[:] = wk1.x.array # Update previous solution with current solution
file.close()