import os
import subprocess
from settingsMaker import Settings

mySettings = Settings()

MaxNbThreads = int(subprocess.check_output(["getconf","_NPROCESSORS_ONLN"]))
print("Max Nb Threads = " + str(MaxNbThreads))

mySettings["finite element solver"]["nbThreads"] = MaxNbThreads
mySettings["demagnetization field solver"]["nbThreads"] = MaxNbThreads

mySettings["outputs"]["evol columns"] = ["t","<mx>","<my>","<mz>","E_tot"]

mySettings["outputs"]["directory"] = "test_data_out/"

mySettings["mesh"]["filename"] = "test.msh"

mySettings["Bext"] = [0.1,0.0,1.0]
mySettings["mesh"]["volume_regions"]["300"]["alpha"] = 1.0

mySettings["time integration"]["final_time"] = 1e-2
mySettings["time integration"]["min(dt)"] = 1e-14

mySettings["initial magnetization"] = {"Mx":"0","My":"0","Mz":"-1"}

mySettings.write('mySettings.json')

os.system("mkdir " + mySettings["outputs"]["directory"])
os.system("./feellgood mySettings.json")

	
