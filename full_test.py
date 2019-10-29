import os
import sys
import subprocess
from math import sqrt
from settingsMaker import Settings

mySettings = Settings()

MaxNbThreads = int(subprocess.check_output(["getconf","_NPROCESSORS_ONLN"]))
print("Max Nb Threads = " + str(MaxNbThreads))

mySettings["finite element solver"]["nbThreads"] = 4 #MaxNbThreads
mySettings["demagnetization field solver"]["nbThreads"] = 4 #MaxNbThreads

mySettings["outputs"]["evol columns"] = ["t","<mx>","<my>","<mz>","E_ex","E_demag","E_zeeman","E_tot"]
mySettings["outputs"]["take_photo"] = 100
mySettings["outputs"]["directory"] = "test_data_out/"

mySettings["mesh"]["filename"] = "ellipsoid.msh"
mySettings["mesh"]["scaling factor"] = 1e-10

mySettings["Bext"] = [1.0,0.0,-1.0]
mySettings["mesh"]["volume_regions"]["300"]["alpha"] = 0.5

mySettings["time integration"]["final_time"] = 1.0e-5
mySettings["time integration"]["min(dt)"] = 1e-8
mySettings["time integration"]["max(dt)"] = 1e-5
mySettings["time integration"]["initial dt"] = 1e-7

mySettings["time integration"]["max(du)"] = 0.2

mySettings["initial magnetization"] = {"Mx":"0","My":"0","Mz":"1"}

mySettings.write('mySettings.json')

if(os.path.exists(mySettings["outputs"]["directory"]) and os.path.isdir(mySettings["outputs"]["directory"]) ):
	print("directory " + mySettings["outputs"]["directory"] + " already exists.")
else:
	os.system("mkdir " + mySettings["outputs"]["directory"])

print("FeeLLGood running...")
val = subprocess.run(["./feellgood","mySettings.json"])

if(val.returncode==0):
	print("FeeLLGood terminated correctly")
	with open(mySettings["outputs"]["directory"] + mySettings["outputs"]["file basename"] + ".evol","r") as f:
		for line in f:
			pass
		lastLine = line
	f.close()
	data = lastLine.split()
	mx = float(data[1])
	my = float(data[2])
	mz = float(data[3])
	print("mag= ",mx,';',my,';',mz)
	X = 0.599620
	Y = -0.000982
	Z = -0.800284
	if(sqrt((X-mx)**2+(Y-my)**2+(Z-mz)**2) < 5e-5):
		valRet = 0
	else:
		valRet = 1
else:
	print("FeeLLGood terminated before final time")
	valRet = 1
print("return " + str(valRet) + " ;dist = ",sqrt((X-mx)**2+(Y-my)**2+(Z-mz)**2))
sys.exit(valRet)

