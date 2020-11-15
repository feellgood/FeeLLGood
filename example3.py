import numpy as np

from icoMaker import Ico
from settingsMaker import Settings

meshFileName = 'icosahedron.msh'
#carefull here, volume and region names must be stringified integers
vol_region_name = "300"
surf_region_name = "200"
my_ico = Ico()
my_ico.make(meshFileName,vol_region_name,surf_region_name)
my_ico.show()

mySettings = Settings(meshFileName)

mySettings["outputs"]["file_basename"] = "ico"
mySettings["outputs"]["evol_time_step"] = 0.2e-8 
mySettings["outputs"]["evol_columns"] = [ "t", "<Mx>", "<My>", "<Mz>", "E_ex", "E_aniso", "E_demag", "E_tot" ]

mySettings.createVolRegion( vol_region_name )
mySettings.createSurfRegion( surf_region_name )

# cubic anisotropy for volume region "300" (Fe : K3 = 4.2e4 J/m^3)
mySettings["mesh"]["volume_regions"]["300"]["K3"] = 4.2e4

#cosine directions [alpha,beta,gamma]
mySettings["mesh"]["volume_regions"]["300"]["alpha"] = [1.0,0.0,0.0]
mySettings["mesh"]["volume_regions"]["300"]["beta"] = [0.0,1.0,0.0]
mySettings["mesh"]["volume_regions"]["300"]["gamma"] = [0.0,0.0,1.0]

mySettings["initial_magnetization"] = {"Mx":"1","My":"0.02","Mz":"0"}

mySettings.write('ico_test.json')

