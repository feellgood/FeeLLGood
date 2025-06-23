import os
import sys

from feellgood.settingsMaker import Settings

mySettings = Settings("ellipsoid.msh")
mySettings.createVolRegion( "ellipsoid_volume" )
mySettings.createSurfRegion( "ellipsoid_surface" )

mySettings["Bext"] = [0.01, 0.0, 0.0]
mySettings["outputs"]["evol_columns"] = ["t","<Mx>","<My>","<Mz>","E_tot"]

mySettings.write('mySettings.json')
sys.stdout.flush()
os.system("../feellgood mySettings.json")
