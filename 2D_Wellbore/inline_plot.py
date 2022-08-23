# This function plots the dolfinx FEM Mesh

import pyvista
pyvista.set_jupyter_backend("pythreejs")
from dolfinx.plot import create_vtk_mesh

def plot_mesh(mesh, cell_tags):
    plotter = pyvista.Plotter()
    grid = pyvista.UnstructuredGrid(*create_vtk_mesh(mesh, mesh.topology.dim))
    num_local_cells = mesh.topology.index_map(mesh.topology.dim).size_local
    grid.cell_data["Marker"] = cell_tags.values[cell_tags.indices<num_local_cells]
    grid.set_active_scalars("Marker")
    actor = plotter.add_mesh(grid, show_edges=True)
    plotter.view_xy()
    if not pyvista.OFF_SCREEN:
        plotter.show()
    else:
        pyvista.start_xvfb()
        cell_tag_fig = plotter.screenshot("cell_tags.png")