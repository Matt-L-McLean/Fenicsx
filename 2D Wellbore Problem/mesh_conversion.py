#-----------------------------------------------------------------------------------
# This script takes a Gmsh API model (python based) as an input and returns the
# Fenicsx equivalent mesh for FEM modeling. This script allows for both 2D and 3D
# models.
#-----------------------------------------------------------------------------------

from dolfinx.graph import create_adjacencylist
from dolfinx.io import (cell_perm_gmsh, distribute_entity_data, extract_gmsh_geometry, 
                        extract_gmsh_topology_and_markers, ufl_mesh_from_gmsh, XDMFFile)
from dolfinx.mesh import CellType, create_mesh, meshtags_from_entities
from mpi4py import MPI
import numpy as np
import gmsh
import warnings

warnings.filterwarnings("ignore")
rank = MPI.COMM_WORLD.rank

def msh_to_xdmf(model, gdim):
    x = extract_gmsh_geometry(model)
    element_types, element_tags, node_tags = model.mesh.getElements(dim=gdim)
    assert len(element_types) == 1
    name, dim, order, num_nodes, local_coords, num_first_order_nodes = model.mesh.getElementProperties(element_types[0])
    cells = node_tags[0].reshape(-1, num_nodes) - 1

    mesh = create_mesh(MPI.COMM_SELF, cells, x, ufl_mesh_from_gmsh(element_types[0], x.shape[1]))
    with XDMFFile(MPI.COMM_SELF, f"Mesh/Gmsh_to_Fenicsx.xdmf", "w") as file:
        file.write_mesh(mesh)
    return mesh