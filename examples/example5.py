import sys
import numpy as np

sys.path.insert(0,'../tools')
from cuboidMaker import Cuboid
from settingsMaker import Settings

rectangle = Cuboid([-2,-2,-2],[2,2,2],4,4,4)
meshFileName = 'rectangle.msh'
vol_region_name = "300"
surf_region_name = "200"

bc1_regNumber = 201
bc2_regNumber = 202

rectangle.add_sub_surface(str(bc1_regNumber), lambda x,y,z: z== -2 )
rectangle.add_sub_surface(str(bc2_regNumber), lambda x,y,z: z== 2 )
rectangle.make(meshFileName,vol_region_name,surf_region_name)

mySettings = Settings(meshFileName)

mySettings["finite_element_solver"]["nb_threads"] = 32 # MaxNbThreads
mySettings["demagnetizing_field_solver"]["nb_threads"] = 32 # MaxNbThreads

mySettings["outputs"]["file_basename"] = "rectangle"

mySettings["outputs"]["evol_time_step"] = 0.5e-11
mySettings["outputs"]["evol_columns"] = [ "t", "<Mx>", "<My>", "<Mz>", "E_ex", "E_aniso", "E_demag", "E_tot" ]

mySettings["outputs"]["take_photo"] = 500

mySettings.createVolRegion( vol_region_name )
mySettings.createSurfRegion( surf_region_name )

#magnetization at saturation (SI unit = A/m)
mySettings["mesh"]["volume_regions"][vol_region_name]["Js"] = 800e3

# exchange constant (unit = J/m)
mySettings["mesh"]["volume_regions"][vol_region_name]["Ae"] = 13e-12

mySettings["mesh"]["volume_regions"][vol_region_name]["alpha_LLG"] = 0.02

mySettings["time_integration"]["final_time"] = 0.5e-9
mySettings["time_integration"]["min(dt)"] = 0.1e-17
mySettings["time_integration"]["max(dt)"] = 1e-9

mySettings["time_integration"]["max(du)"] = 0.001

gamma0 = 1
sigma = 1
N0 = 1
beta = 1
l_J = 1
l_sf = 1
V1 = 1.2345
V2 = -3.14
mySettings.createSTT(300,gamma0,sigma,N0,beta,l_J,l_sf,bc1_regNumber,"V",V1,bc2_regNumber,"V",V2)

mySettings["initial_magnetization"] = {"Mx":"x","My":"y","Mz":"0.1"}

mySettings.write('rectangleSTT_test.json')
