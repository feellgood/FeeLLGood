import os
from settingsMaker import Settings

mySettings = Settings()

print("output directory is " + mySettings["outputs"]["directory"] )

myFields = [0.01,0.02,0.04]

for Bz in myFields:
	myDir = "data_out_B" + str(Bz) + "/"
	mySettings["outputs"]["directory"] = myDir
	os.system("mkdir "+ myDir)
	mySettings["Bext"] = [0.0,0.0,Bz]	
	print("Bext is " + str(toto["Bext"]))
	mySettings.write('mySettings.json')
	os.system("./feellgood mySettings.json")
