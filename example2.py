import os
import sys
import subprocess
from math import sqrt
from settingsMaker import Settings

mySettings = Settings()

MaxNbThreads = int(subprocess.check_output(["getconf","_NPROCESSORS_ONLN"]))

mySettings["finite_element_solver"]["nb_threads"] = 4 # MaxNbThreads
mySettings["finite_element_solver"]["max(iter)"] = 700

mySettings["demagnetization_field_solver"]["nb_threads"] = 4 # MaxNbThreads

mySettings["outputs"]["evol_columns"] = ["t","<Mx>","<My>","<Mz>","<dMx/dt>","<dMy/dt>","<dMz/dt>","Hx","Hy","Hz","E_ex","E_demag","E_zeeman","E_tot"]
mySettings["outputs"]["take_photo"] = 500
mySettings["outputs"]["directory"] = "test_data_out/"
mySettings["outputs"]["verbose"] = False

mySettings["outputs"]["evol_time_step"] = 0.1e-9

mySettings["mesh"]["filename"] = "ellipsoid.msh"
mySettings["mesh"]["scaling_factor"] = 1e-10

A = 0.01
f = 5e8

mySettings["Bext"] = {"Bx" : str(A) + "*cos(2*Pi*" + str(f) + "*t)", "By" : str(A) + "*sin(2*Pi*"  + str(f) +  "*t)" , "Bz": "10"}
mySettings["mesh"]["volume_regions"]["300"]["alpha"] = 0.05

mySettings["time_integration"]["final_time"] = 5.0e-7
mySettings["time_integration"]["min(dt)"] = 1e-12
mySettings["time_integration"]["max(dt)"] = 5e-9
mySettings["time_integration"]["initial_dt"] = 0.5e-10

mySettings["time_integration"]["max(du)"] = 0.1

mySettings["initial_magnetization"] = {"Mx":"1","My":"1","Mz":"1"}

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
	print("m_rho = %.2e"% sqrt((mx)**2+(my)**2))
sys.exit(0)

