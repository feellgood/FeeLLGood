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

mySettings.write('mySettings.json')

os.system("mkdir " + mySettings["outputs"]["directory"])
os.system("./feellgood mySettings.json")

	
