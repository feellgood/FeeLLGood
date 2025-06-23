import os
import sys

from feellgood.settingsMaker import Settings

mySettings = Settings("ellipsoid.msh")
mySettings.createVolRegion( "ellipsoid_volume" )
mySettings.createSurfRegion( "ellipsoid_surface" )

myFields = {"B0" : 0.01, "B1" : 0.02 , "B2": 0.04}

for B in myFields:
    myDir = "data_out_B" + str(myFields[B]) + "/"
    mySettings["outputs"]["directory"] = myDir
    os.system("mkdir "+ myDir)
    mySettings["Bext"] = [0.0, 0.0, myFields[B]]
    print("Bext is " + str(mySettings["Bext"]))
    mySettings.write('mySettings.json')
    sys.stdout.flush()
    os.system("../feellgood mySettings.json")
