import sys
import numpy as np

from cuboidMaker import Cuboid

petitCube = Cuboid([-2,-2,-2],[2,2,2],2,2,2)

meshFileName = 'petitCube.msh'
vol_region_name = "300"
surf_region_name = "200"

#func = lambda x,y,z: z== -2
petitCube.add_sub_surface("201", lambda x,y,z: z== -2 )

petitCube.add_sub_surface("202", lambda x,y,z: z== 2 )

petitCube.make(meshFileName,vol_region_name,surf_region_name)

