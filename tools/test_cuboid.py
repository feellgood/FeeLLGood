import sys
import numpy as np

from cuboidMaker import Cuboid

zmin = -3
zmax = 3

petitCube = Cuboid([-2,-2,zmin],[2,2,zmax],20,20,20)

meshFileName = 'petitCube.msh'
vol_region_name = "300"
surf_region_name = "200"

#func = lambda x,y,z: z== -2
petitCube.add_sub_surface("201", lambda x,y,z: z== zmin )

petitCube.add_sub_surface("202", lambda x,y,z: z== zmax )

# z0 is the plane of symetry of the modulation transverse to the axis of the wire
z0 = 1.0

# A is the amplitude of the modulation, if negative it is a constriction, if positive a protrusion
A = -0.2
def fx(x,y,z):
    return x*(1+A/(1+(z-z0)**2))

def fy(x,y,z):
    return y*(1+A/(1+(z-z0)**2))

def fz(x,y,z):
    return z

petitCube.transformNodes(lambda pt:[fx(pt[0],pt[1],pt[2]),fy(pt[0],pt[1],pt[2]),fz(pt[0],pt[1],pt[2])])

petitCube.make(meshFileName,vol_region_name,surf_region_name)

