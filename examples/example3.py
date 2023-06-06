import sys
import numpy as np

from feellgood.meshMaker import Ico
from feellgood.settingsMaker import Settings

meshFileName = 'icosahedron.msh'
vol_region_name = "whole_volume"
surf_region_name = "whole_surface"
my_ico = Ico()
my_ico.make(meshFileName,vol_region_name,surf_region_name)

mySettings = Settings(meshFileName)

mySettings["outputs"]["file_basename"] = "ico"

mySettings["outputs"]["evol_time_step"] = 1e-17
mySettings["outputs"]["final_time"] = 5e-14
mySettings["outputs"]["evol_columns"] = [ "t", "<Mx>", "<My>", "<Mz>", "E_ex", "E_aniso", "E_demag","E_zeeman", "E_tot" ]

mySettings["outputs"]["take_photo"] = 1000

mySettings.createVolRegion( vol_region_name )
mySettings.createSurfRegion( surf_region_name )

# cubic anisotropy for volume region "whole_volume" (Fe : K3 = 4.2e4 J/m^3)
mySettings["mesh"]["volume_regions"]["whole_volume"]["K3"] = 4.2e4

#magnetization at saturation (SI unit = A/m)
mySettings["mesh"]["volume_regions"]["whole_volume"]["Js"] = 800e3

# exchange constant (unit = J/m)
mySettings["mesh"]["volume_regions"]["whole_volume"]["Ae"] = 10e-12

mySettings["mesh"]["volume_regions"]["whole_volume"]["alpha_LLG"] = 0.05

mySettings["time_integration"]["min(dt)"] = 5e-18
mySettings["time_integration"]["max(dt)"] = 5e-16

mySettings["time_integration"]["max(du)"] = 40


#cosine directions [alpha,beta,gamma]
mySettings["mesh"]["volume_regions"]["whole_volume"]["ex"] = [1.0,0.0,0.0]
mySettings["mesh"]["volume_regions"]["whole_volume"]["ey"] = [0.0,1.0,0.0]
mySettings["mesh"]["volume_regions"]["whole_volume"]["ez"] = [0.0,0.0,1.0]

mySettings["initial_magnetization"] = [1, 0, 0]

mySettings.write('ico_test.json')

