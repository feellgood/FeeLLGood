import os
from settingsMaker import Settings

mySettings = Settings("ellipsoid.msh")
mySettings.createVolRegion( "300" )
mySettings.createSurfRegion( "200" )

print("output directory is " + mySettings["outputs"]["directory"] )

myFields = {"B0" : "0.01", "B1" : "0.02" , "B2": "0.04"}

for B in myFields:
	myDir = "data_out_B" + str(myFields[B]) + "/"
	mySettings["outputs"]["directory"] = myDir
	os.system("mkdir "+ myDir)
	mySettings["Bext"] = {"Bx" : "0.0", "By" : "0.0" , "Bz": str(myFields[B])}	
	print("Bext is " + str(mySettings["Bext"]))
	mySettings.write('mySettings.json')
	os.system("./feellgood mySettings.json")
