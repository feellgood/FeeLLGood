import os
from settingsMaker import Settings

toto = Settings()

print("output directory is " + toto["outputs"]["directory"] )

myFields = [0.01,0.02,0.04]

for Bz in myFields:
	myDir = "data_out_B" + str(Bz) + "/"
	toto["outputs"]["directory"] = myDir
	os.system("mkdir "+ myDir)
	toto["Bext"] = [0.0,0.0,Bz]	
	print("Bext is " + str(toto["Bext"]))
	toto.write('essaiSettings.json')
	os.system("./feellgood essaiSettings.json")
