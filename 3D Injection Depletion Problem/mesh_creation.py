#-----------------------------------------------------------------------------------
# This script takes model dimensions and mesh sizes and returns the Gmsh API model
# and model dimension, e.g., 2D or 3D. The mesh employs three-dimensional tetrahedral
# elements for volumes and triangles for surfaces.
#-----------------------------------------------------------------------------------

import gmsh
from mpi4py import MPI

rank = MPI.COMM_WORLD.rank

def Gmsh_model(xl, yl, zl, ms, msw):
    gmsh.model.setCurrent("Reservoir")
    gmsh.initialize()
    gdim = 3
    if rank == 0:
        #-----------------------------------------
        # boundary nodes
        #-----------------------------------------
        p1 = gmsh.model.geo.addPoint(0,0,0,ms,1)
        p2 = gmsh.model.geo.addPoint(xl,0,0,ms,2)
        p3 = gmsh.model.geo.addPoint(xl,yl,0,ms,3)
        p4 = gmsh.model.geo.addPoint(0,yl,0,ms,4)
        p5 = gmsh.model.geo.addPoint(0,0,zl,ms,5)
        p6 = gmsh.model.geo.addPoint(xl,0,zl,ms,6)
        p7 = gmsh.model.geo.addPoint(xl,yl,zl,ms,7)
        p8 = gmsh.model.geo.addPoint(0,yl,zl,ms,8)
        #-----------------------------------------
        # boundary lines
        #-----------------------------------------
        l1 = gmsh.model.geo.addLine(p1,p2,1)
        l2 = gmsh.model.geo.addLine(p2,p3,2)
        l3 = gmsh.model.geo.addLine(p3,p4,3)
        l4 = gmsh.model.geo.addLine(p4,p1,4)
        l5 = gmsh.model.geo.addLine(p5,p6,5)
        l6 = gmsh.model.geo.addLine(p6,p7,6)
        l7 = gmsh.model.geo.addLine(p7,p8,7)
        l8 = gmsh.model.geo.addLine(p8,p5,8)
        l9 = gmsh.model.geo.addLine(p1,p5,9)
        l10 = gmsh.model.geo.addLine(p2,p6,10)
        l11 = gmsh.model.geo.addLine(p3,p7,11)
        l12 = gmsh.model.geo.addLine(p4,p8,12)
        #-----------------------------------------
        # boundary curve loop
        #-----------------------------------------
        cl1 = gmsh.model.geo.addCurveLoop([l1,l10,-l5,-l9],1)
        cl2 = gmsh.model.geo.addCurveLoop([l2,l11,-l6,-l10],2)
        cl3 = gmsh.model.geo.addCurveLoop([l3,l12,-l7,-l11],3)
        cl4 = gmsh.model.geo.addCurveLoop([l4,l9,-l8,-l12],4)
        cl5 = gmsh.model.geo.addCurveLoop([l1,l2,l3,l4],5)
        cl6 = gmsh.model.geo.addCurveLoop([l5,l6,l7,l8],6)
        #-----------------------------------------
        # boundary surface
        #-----------------------------------------
        s1 = gmsh.model.geo.addPlaneSurface([cl1])
        s2 = gmsh.model.geo.addPlaneSurface([cl2])
        s3 = gmsh.model.geo.addPlaneSurface([cl3])
        s4 = gmsh.model.geo.addPlaneSurface([cl4])
        s5 = gmsh.model.geo.addPlaneSurface([cl5])
        s6 = gmsh.model.geo.addPlaneSurface([cl6])
        #-----------------------------------------
        # surface loop --> combine boundary surfaces
        #-----------------------------------------
        sl1 = gmsh.model.geo.addSurfaceLoop([s1,s2,s3,s4,s5,s6])
        #-----------------------------------------
        # volume
        #-----------------------------------------
        v1 = gmsh.model.geo.addVolume([sl1])
        #-----------------------------------------
        # embedded injection nodes
        #-----------------------------------------
        p9 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9,msw,9)
        p10 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*1,msw,10)
        p11 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*2,msw,11)
        p12 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*3,msw,12)
        p13 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*4,msw,13)
        p14 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*5,msw,14)
        p15 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*6,msw,15)
        p16 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*7,msw,16)
        p17 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*8,msw,17)
        p18 = gmsh.model.geo.addPoint(xl/2,yl/2,(zl-10)/9+80/9*9,msw,18)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(0,[p9,p10,p11,p12,p13,p14,p15,p16,p17,p18],3,v1)
        #-----------------------------------------
        # mesh refinement around wellbore
        #-----------------------------------------
        well = gmsh.model.mesh.field.add("Box",1)
        gmsh.model.mesh.field.setNumber(1, "VIn", msw)
        gmsh.model.mesh.field.setNumber(1, "VOut", ms)
        gmsh.model.mesh.field.setNumber(1, "XMin", xl/2-25)
        gmsh.model.mesh.field.setNumber(1, "XMax", xl/2+25)
        gmsh.model.mesh.field.setNumber(1, "YMin", yl/2-25)
        gmsh.model.mesh.field.setNumber(1, "YMax", yl/2+25)
        gmsh.model.mesh.field.setNumber(1, "ZMin", 0)
        gmsh.model.mesh.field.setNumber(1, "ZMax", zl)
        gmsh.model.mesh.field.setNumber(1, "Thickness", 100)
        gmsh.model.mesh.field.setAsBackgroundMesh(1)
        #-----------------------------------------
        # synchronize
        #-----------------------------------------
        gmsh.model.mesh.generate(gdim)
        return gmsh.model, gdim
