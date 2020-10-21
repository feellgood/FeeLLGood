import os
import sys
import subprocess
from math import sqrt
from settingsMaker import Settings

mySettings = Settings("ellipsoid.msh")

MaxNbThreads = int(subprocess.check_output(["getconf","_NPROCESSORS_ONLN"]))

mySettings["finite_element_solver"]["nb_threads"] = MaxNbThreads
mySettings["finite_element_solver"]["max(iter)"] = 500

mySettings["demagnetization_field_solver"]["nb_threads"] = MaxNbThreads

mySettings["outputs"]["evol_columns"] = ["t","<Mx>","<My>","<Mz>","E_ex","E_demag","E_zeeman","E_tot"]
mySettings["outputs"]["take_photo"] = 100
mySettings["outputs"]["directory"] = "test_data_out/"

mySettings["outputs"]["evol_time_step"] = 1e-7

mySettings["mesh"]["scaling_factor"] = 1e-10

mySettings["Bext"] = {"Bx" : "1", "By" : "0" , "Bz": "-1"}

#carefull here, volume and region names must be stringified integers
vol_region_name = "300"
surf_region_name = "200"

mySettings.createVolRegion( vol_region_name )
mySettings.createSurfRegion( surf_region_name )
mySettings["mesh"]["volume_regions"][vol_region_name]["alpha"] = 0.5

mySettings["time_integration"]["final_time"] = 1.5e-5
mySettings["time_integration"]["min(dt)"] = 1e-11
mySettings["time_integration"]["max(dt)"] = 0.5e-7
mySettings["time_integration"]["initial_dt"] = 0.5e-8

mySettings["time_integration"]["max(du)"] = 0.1

mySettings["initial_magnetization"] = {"Mx":"0","My":"0","Mz":"1"}

mySettings.write('mySettings.json')

if(os.path.exists(mySettings["outputs"]["directory"]) and os.path.isdir(mySettings["outputs"]["directory"]) ):
	print("directory " + mySettings["outputs"]["directory"] + " already exists.")
else:
	os.system("mkdir " + mySettings["outputs"]["directory"])

val = subprocess.run(["./feellgood","mySettings.json"])

if(val.returncode==0):
	print("FeeLLGood terminated correctly")
	with open(mySettings["outputs"]["directory"] + mySettings["outputs"]["file_basename"] + ".evol","r") as f:
		for line in f:
			pass
		lastLine = line
	f.close()
	data = lastLine.split()
	mx = float(data[1])
	my = float(data[2])
	mz = float(data[3])
	print("mag= ",mx,';',my,';',mz)
	X = 0.600566
	Y = -0.00193843
	Z = -0.799573
	if(sqrt((X-mx)**2+(Y-my)**2+(Z-mz)**2) < 1e-5):
		valRet = 0
	else:
		valRet = 1
else:
	print("FeeLLGood terminated before final time")
	valRet = 1
print("return " + str(valRet) + " ;dist = %.2e"% sqrt((X-mx)**2+(Y-my)**2+(Z-mz)**2))
sys.exit(valRet)

