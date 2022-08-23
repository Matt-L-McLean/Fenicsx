# This script takes a Gmsh .msh file as an input and
# returns the dolfinx equivalent mesh for FEM modeling
# This script allows for 2D and 3D meshing

from dolfinx.io import (cell_perm_gmsh, distribute_entity_data, extract_gmsh_geometry, 
                        extract_gmsh_topology_and_markers, ufl_mesh_from_gmsh)
from dolfinx.cpp.mesh import to_type
from dolfinx.graph import create_adjacencylist
from dolfinx.mesh import create_mesh, meshtags_from_entities
from mpi4py import MPI
import numpy as np
import gmsh
import warnings

warnings.filterwarnings("ignore")
rank = MPI.COMM_WORLD.rank

def msh_to_xdmf(model, gdim):
    if rank == 0:
        x = extract_gmsh_geometry(model)
        topologies = extract_gmsh_topology_and_markers(model)
        num_cell_types = len(topologies.keys())
        cell_information = {}
        cell_dimensions = np.zeros(num_cell_types, dtype=np.int32)
        for i, element in enumerate(topologies.keys()):
            properties = model.mesh.getElementProperties(element)
            name, dim, order, num_nodes, local_coords, _ = properties
            cell_information[i] = {"id": element, "dim": dim, "num_nodes": num_nodes}
            cell_dimensions[i] = dim
        # Sort elements by ascending dimension
        perm_sort = np.argsort(cell_dimensions)

        # Broadcast cell type data and geometric dimension
        cell_id = cell_information[perm_sort[-1]]["id"]
        tdim = cell_information[perm_sort[-1]]["dim"]
        num_nodes = cell_information[perm_sort[-1]]["num_nodes"]
        cell_id, num_nodes = MPI.COMM_WORLD.bcast([cell_id, num_nodes], root=0)

        cells = np.asarray(topologies[cell_id]["topology"], dtype=np.int64)
        cell_values = np.asarray(topologies[cell_id]["cell_data"], dtype=np.int32)
    else:
        cell_id, num_nodes = MPI.COMM_WORLD.bcast([None, None], root=0)
        cells, x = np.empty([0, num_nodes], dtype=np.int64), np.empty([0, gdim])
        cell_values = np.empty((0,), dtype=np.int32)
    gmsh.finalize()

    # Create distributed mesh
    ufl_domain = ufl_mesh_from_gmsh(cell_id, gdim)
    gmsh_cell_perm = cell_perm_gmsh(to_type(str(ufl_domain.ufl_cell())), num_nodes)
    cells = cells[:, gmsh_cell_perm]
    mesh = create_mesh(MPI.COMM_WORLD, cells, x[:, :gdim], ufl_domain)
    tdim = mesh.topology.dim

    local_entities, local_values = distribute_entity_data(mesh, tdim, cells, cell_values)
    mesh.topology.create_connectivity(tdim, 0)
    adj = create_adjacencylist(local_entities)
    ct = meshtags_from_entities(mesh, tdim, adj, np.int32(local_values))
    return mesh, ct

