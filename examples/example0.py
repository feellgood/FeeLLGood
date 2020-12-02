import os
import sys

sys.path.insert(0,'../tools')
from settingsMaker import Settings

mySettings = Settings("../ellipsoid.msh")
mySettings.createVolRegion( "300" )
mySettings.createSurfRegion( "200" )

mySettings["Bext"] = {"Bx" : "0.01", "By" : "0.0" , "Bz": "0.0"}
mySettings["outputs"]["evol_columns"] = ["t","<Mx>","<My>","<Mz>","E_tot"]

mySettings.write('mySettings.json')
os.system("../feellgood mySettings.json")
