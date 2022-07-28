# Fenicsx
This repository contains 2D and 3D poroelastic numerical simulations solved with open source software Fenicsx/Dolfinx. FEM mesh is generated through gmsh software (alternatively pygmsh with python) and converted to xdmf file for use in Fenicsx. A mixed finite element function space is employed to solve for displacement $\vec{u}=\left(u,v,w\right)$ and pore fluid pressure $p$. Vector piecewise quadratic Lagrange elements are used for displacement and scalar piecewise linear Lagrange elements are used for fluid pressure. The backward Euler method of time integration is used to discretize pore fluid diffusion.

The governing equations of solid equilibrium and mass balance are (Cheng, 2016): $$G\nabla^2\vec{u}+\frac{G}{1-2\nu}\nabla\left(\nabla\cdot\vec{u}\right)-\alpha\nabla{p}=0$$ $$\frac{1}{M^*}\frac{\partial p}{\partial t} - \kappa \nabla^2p + \alpha\frac{\partial \left(\nabla\cdot\vec{u}\right)}{\partial t}=0$$

where $G$ is the shear modulus, $\nu$ is the Poisson's ratio, $\alpha$ is the Biot coefficient, $M^*$ is the Biot modulus, $\kappa$ is the pore fluid mobility (defined as ratio of permeability to viscosity).

Note: Fenicsx is currently unable to handle mixed-dimensional models with discontinuous "jumps" in the displacement field: $[[\vec{u}]]=\vec{u}^+-\vec{u}^-$. Hence, this repository is unable to simulate discrete fractures modeled in $n-1$ dimensions (as is typically done in hydraulic fracturing, geothermal fracture newtorks, etc.).

Refs:
[1] Cheng, A. H. D. (2016). Poroelasticity (Vol. 27). Switzerland: Springer International Publishing.
