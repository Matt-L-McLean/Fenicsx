# This script takes inputs for domain half-length
# in x and y directions, well diameter, and local
# mesh size. This script returns the gmsh.model file

import gmsh
from mpi4py import MPI

rank = MPI.COMM_WORLD.rank

def Gmsh_model(ms, mw, xc, yc, wc):
    gmsh.initialize()
    gdim = 2
    if rank == 0:

        # boundary nodes
        ll = gmsh.model.geo.addPoint(-xc,-yc,0,ms,1)
        lr = gmsh.model.geo.addPoint(xc,-yc,0,ms,2)
        ur = gmsh.model.geo.addPoint(xc,yc,0,ms,3)
        ul = gmsh.model.geo.addPoint(-xc,yc,0,ms,4)
        cen = gmsh.model.geo.addPoint(0,0,0,mw,5)
        uw = gmsh.model.geo.addPoint(0,wc/2,0,mw,6)
        lw = gmsh.model.geo.addPoint(0,-wc/2,0,mw,7)

        # boundary lines
        bt = gmsh.model.geo.addLine(ll,lr,1)
        rt = gmsh.model.geo.addLine(lr,ur,2)
        tp = gmsh.model.geo.addLine(ur,ul,3)
        lt = gmsh.model.geo.addLine(ul,ll,4)
        rw = gmsh.model.geo.addCircleArc(lw,cen,uw,5)
        rl = gmsh.model.geo.addCircleArc(uw,cen,lw,6)

        # boundary curve loop
        bd = gmsh.model.geo.addCurveLoop([bt,rt,tp,lt],1)
        ib = gmsh.model.geo.addCurveLoop([rw,rl],2)

        # boundary surface
        s = gmsh.model.geo.addPlaneSurface([bd,ib])

        # physical tags for BCs
        base = gmsh.model.addPhysicalGroup(1, [bt])
        gmsh.model.setPhysicalName(1, base, "base")

        right = gmsh.model.addPhysicalGroup(1, [rt])
        gmsh.model.setPhysicalName(1, right, "right")

        left = gmsh.model.addPhysicalGroup(1, [lt])
        gmsh.model.setPhysicalName(1, left, "left")

        top = gmsh.model.addPhysicalGroup(1, [tp])
        gmsh.model.setPhysicalName(1, top, "top")

        well = gmsh.model.addPhysicalGroup(1, [ib])
        gmsh.model.setPhysicalName(1, well, "well")

        surf = gmsh.model.addPhysicalGroup(2, [s])
        gmsh.model.setPhysicalName(2, s, "surf")

        # synchronize
        gmsh.model.geo.synchronize()
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(gdim)
        return gmsh.model, gdim