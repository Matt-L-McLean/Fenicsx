# Poroelasticity: Numerical Simulation
This repository contains 2D and 3D poroelastic numerical simulations solved with open source software Fenicsx/Dolfinx (Scroggs et al., 2022; Alnaes et al., 2014). FEM mesh is generated through gmsh API software (python) and converted to .xdmf file for use in Fenicsx. Post-processing/results visualization is with Paraview software.

The governing equations of solid equilibrium and mass balance are (Cheng, 2016): 

<img src="https://latex.codecogs.com/svg.image?G\nabla^2\vec{u}&plus;\frac{G}{1-2\nu}\nabla\left(\nabla\cdot\vec{u}\right)-\alpha\nabla{p}=0" title="https://latex.codecogs.com/svg.image?G\nabla^2\vec{u}+\frac{G}{1-2\nu}\nabla\left(\nabla\cdot\vec{u}\right)-\alpha\nabla{p}=0" />

<img src="https://latex.codecogs.com/svg.image?\frac{1}{M^*}\frac{\partial&space;p}{\partial&space;t}&space;-&space;\kappa&space;\nabla^2p&space;&plus;&space;\alpha\frac{\partial&space;\left(\nabla\cdot\vec{u}\right)}{\partial&space;t}=0" title="https://latex.codecogs.com/svg.image?\frac{1}{M^*}\frac{\partial p}{\partial t} - \kappa \nabla^2p + \alpha\frac{\partial \left(\nabla\cdot\vec{u}\right)}{\partial t}=0" />

where <img src="https://latex.codecogs.com/svg.image?\inline&space;G" title="https://latex.codecogs.com/svg.image?\inline G" /> is the shear modulus, <img src="https://latex.codecogs.com/svg.image?\inline&space;\nu" title="https://latex.codecogs.com/svg.image?\inline \nu" /> is the Poisson's ratio, <img src="https://latex.codecogs.com/svg.image?\inline&space;\alpha" title="https://latex.codecogs.com/svg.image?\inline \alpha" /> is the Biot coefficient, <img src="https://latex.codecogs.com/svg.image?\inline&space;M^*" title="https://latex.codecogs.com/svg.image?\inline M^*" /> is the Biot modulus, <img src="https://latex.codecogs.com/svg.image?\inline&space;\kappa=k/\mu" title="https://latex.codecogs.com/svg.image?\inline \kappa=k/\mu" /> is the pore fluid mobility (defined as ratio of permeability to viscosity).

Note: Fenicsx is currently unable to handle mixed-dimensional models with discontinuous "jumps" in the displacement field: <img src="https://latex.codecogs.com/svg.image?\inline&space;[[\vec{u}]]=\vec{u}^&plus;-\vec{u}^-" title="https://latex.codecogs.com/svg.image?\inline [[\vec{u}]]=\vec{u}^+-\vec{u}^-" />. Hence, this repository is unable to simulate discrete fractures modeled in <img src="https://latex.codecogs.com/svg.image?\inline&space;n-1" title="https://latex.codecogs.com/svg.image?\inline n-1" /> dimensions (as is typically done in hydraulic fracturing, geothermal fracture newtorks, etc.).

Refs:

[1] M. W. Scroggs, J. S. Dokken, C. N. Richardson, and G. N. Wells. Construction of arbitrary order finite element degree-of-freedom maps on polygonal and polyhedral cell meshes, ACM Transactions on Mathematical Software (2022).

[2] M. S. Alnaes, A. Logg, K. B. Ølgaard, M. E. Rognes and G. N. Wells. Unified Form Language: A domain-specific language for weak formulations of partial differential equations, ACM Transactions on Mathematical Software 40 (2014).

[3] Cheng, A. H. D. Poroelasticity (Vol. 27). Switzerland: Springer International Publishing (2016).
