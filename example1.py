import os
from settingsMaker import Settings

mySettings = Settings()

print("output directory is " + mySettings["outputs"]["directory"] )

myFields = {"Bx" : "0.01", "By" : "0.02" , "Bz": "0.04"}

for Bz in myFields:
	myDir = "data_out_B" + str(Bz) + "/"
	mySettings["outputs"]["directory"] = myDir
	os.system("mkdir "+ myDir)
	mySettings["Bext"] = {"Bx" : "0.0", "By" : "0.0" , "Bz": str(Bz)}	
	print("Bext is " + str(mySettings["Bext"]))
	mySettings.write('mySettings.json')
	os.system("./feellgood mySettings.json")
