# Poroelasticity Simulation
This repository contains 2D and 3D poroelastic numerical simulations solved with open source software Fenicsx/Dolfinx (Scroggs et al., 2022; Alnaes et al., 2014). FEM mesh is generated through gmsh software (alternatively pygmsh with python) and converted to xdmf file for use in Fenicsx. A mixed finite element function space is employed to solve for displacement $\vec{u}=\left(u,v,w\right)$ and pore fluid pressure $p$. Vector piecewise quadratic Lagrange elements are used for displacement and scalar piecewise linear Lagrange elements are used for fluid pressure. The Euler method of time integration is used to discretize pore fluid diffusion.

A geomechanics sign covention is employed: compressive stress and contraction strain are positive. Units are pascal for stress and meter for length.

The governing equations of solid equilibrium and mass balance are (Cheng, 2016): $$G\nabla^2\vec{u}+\frac{G}{1-2\nu}\nabla\left(\nabla\cdot\vec{u}\right)-\alpha\nabla{p}=0$$ $$\frac{1}{M^*}\frac{\partial p}{\partial t} - \kappa \nabla^2p + \alpha\frac{\partial \left(\nabla\cdot\vec{u}\right)}{\partial t}=0$$

where $G$ is the shear modulus, $\nu$ is the Poisson's ratio, $\alpha$ is the Biot coefficient, $M^*$ is the Biot modulus, $\kappa$ is the pore fluid mobility (defined as ratio of permeability to viscosity).

Note: Fenicsx is currently unable to handle mixed-dimensional models with discontinuous "jumps" in the displacement field: $[[\vec{u}]]=\vec{u}^+-\vec{u}^-$. Hence, this repository is unable to simulate discrete fractures modeled in $n-1$ dimensions (as is typically done in hydraulic fracturing, geothermal fracture newtorks, etc.).

The poroelastic equations in the weak form with variational test functions $\delta u$ and $\delta p$ are: $$-\int_{\Omega}\nabla\cdot\left(\sigma^{k+1}:\delta \varepsilon^{k+1}\right)\mathrm{d}\Omega=0$$ $$\int_{\Omega}\frac{\alpha^2}{K_u-K_d}\left(\frac{p^{k+1}-p^k}{\Delta t}\right)\delta p^{k+1}\mathrm{d}\Omega - \int_{\Omega}\kappa\left(\nabla p^{k+1}\cdot \nabla\delta p^{k+1} \right)\mathrm{d}\Omega + \int_{\Omega}\frac{\alpha}{\Delta t}\left(\nabla\cdot u^{k+1}-\nabla\cdot u^{k}\right)\delta p^{k+1}\mathrm{d}\Omega=0$$

where $\sigma$ is the total stress, $\varepsilon = 0.5\left(\nabla u + \nabla^T u\right)$ is strain, $\Delta t$ is the time increment, $K_u$ is the undrained bulk modulus, and $K_d$ is the drained bulk modulus. Displacements and fluid pressure with $k+1$ superscript are the current solution while those with $k$ superscript are the solution from previous time step. Other methods of time discretization could be used for increased solution accuracy (e.g., Crank-Nicolson, backward Euler, etc.), however the Euler method is simple to implement.

Refs:

[1] M. W. Scroggs, J. S. Dokken, C. N. Richardson, and G. N. Wells. Construction of arbitrary order finite element degree-of-freedom maps on polygonal and polyhedral cell meshes, ACM Transactions on Mathematical Software (2022).

[2] M. S. Alnaes, A. Logg, K. B. Ølgaard, M. E. Rognes and G. N. Wells. Unified Form Language: A domain-specific language for weak formulations of partial differential equations, ACM Transactions on Mathematical Software 40 (2014).

[3] Cheng, A. H. D. Poroelasticity (Vol. 27). Switzerland: Springer International Publishing (2016).
