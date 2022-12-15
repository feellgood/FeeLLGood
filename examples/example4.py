import sys
import numpy as np

sys.path.insert(0,'../tools')
from cuboidMaker import Cuboid
from settingsMaker import Settings

rectangle = Cuboid([-256,-256,-2],[256,256,2],128,128,1)
meshFileName = 'rectangle.msh'
vol_region_name = "300"
surf_region_name = "200"
rectangle.make(meshFileName,vol_region_name,surf_region_name)

mySettings = Settings(meshFileName)

mySettings["finite_element_solver"]["nb_threads"] = 32 # MaxNbThreads
mySettings["demagnetizing_field_solver"]["nb_threads"] = 32 # MaxNbThreads

mySettings["outputs"]["file_basename"] = "rectangle"

mySettings["outputs"]["evol_time_step"] = 2.5e-17
mySettings["outputs"]["evol_columns"] = [ "t", "<Mx>", "<My>", "<Mz>", "E_ex", "E_aniso", "E_demag", "E_tot" ]

mySettings["outputs"]["take_photo"] = 500

mySettings.createVolRegion( vol_region_name )
mySettings.createSurfRegion( surf_region_name )

# cubic anisotropy for volume region "300" (Fe : K3 = 4.2e4 J/m^3)
mySettings["mesh"]["volume_regions"]["300"]["K3"] = -4.2e4

#magnetization at saturation (SI unit = A/m)
mySettings["mesh"]["volume_regions"]["300"]["Js"] = 800e3

# exchange constant (unit = J/m)
mySettings["mesh"]["volume_regions"]["300"]["Ae"] = 13e-12

mySettings["mesh"]["volume_regions"]["300"]["alpha_LLG"] = 0.02

mySettings["time_integration"]["final_time"] = 2.5e-15
mySettings["time_integration"]["min(dt)"] = 5e-24
mySettings["time_integration"]["max(dt)"] = 5e-15

mySettings["time_integration"]["max(du)"] = 0.001


#cosine directions [alpha,beta,gamma]
mySettings["mesh"]["volume_regions"]["300"]["ex"] = [1.0,0.0,0.0]
mySettings["mesh"]["volume_regions"]["300"]["ey"] = [0.0,1.0,0.0]
mySettings["mesh"]["volume_regions"]["300"]["ez"] = [0.0,0.0,1.0]

mySettings["initial_magnetization"] = ["x", "y", "0.1"]

mySettings.write('rectangle_test.json')
