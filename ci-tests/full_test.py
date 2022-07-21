#!/usr/bin/python3

import os
import sys
import subprocess
from math import sqrt

# Move to the script's directory.
os.chdir(sys.path[0])

sys.path.insert(0,'../tools')
from settingsMaker import Settings

mySettings = Settings("../examples/ellipsoid.msh")

MaxNbThreads = int(subprocess.check_output(["getconf","_NPROCESSORS_ONLN"]))
nbThreads = (MaxNbThreads + 4) // 3 # ellipsoid.msh is small, maxnbthreads is counterproductive

mySettings["finite_element_solver"]["nb_threads"] = nbThreads
mySettings["finite_element_solver"]["max(iter)"] = 500

mySettings["demagnetizing_field_solver"]["nb_threads"] = nbThreads

mySettings["outputs"]["evol_columns"] = ["t","<Mx>","<My>","<Mz>","E_ex","E_demag","E_zeeman","E_tot"]
mySettings["outputs"]["take_photo"] = 100
mySettings["outputs"]["directory"] = "test_data_out/"
mySettings["outputs"]["file_basename"] = "ellipsoid"

mySettings["outputs"]["evol_time_step"] = 1e-7

mySettings["mesh"]["scaling_factor"] = 1e-10

mySettings["Bext"] = [1, 0, -1]

#carefull here, volume and region names must be stringified integers
vol_region_name = "300"
surf_region_name = "200"

mySettings.createVolRegion( vol_region_name )
mySettings.createSurfRegion( surf_region_name )
mySettings["mesh"]["volume_regions"][vol_region_name]["alpha_LLG"] = 0.5

mySettings["time_integration"]["final_time"] = 1.5e-5
mySettings["time_integration"]["min(dt)"] = 1e-11
mySettings["time_integration"]["max(dt)"] = 0.5e-7

mySettings["time_integration"]["max(du)"] = 0.1

mySettings["initial_magnetization"] = [0, 0, 1]

JSON_fileName = 'full_test_settings.json'
mySettings.write(JSON_fileName)

val = subprocess.run(["../feellgood",JSON_fileName])

if(val.returncode==0):
	with open(mySettings["outputs"]["directory"] + mySettings["outputs"]["file_basename"] + ".evol","r") as f:
		for line in f:
			pass
		lastLine = line
	f.close()
	data = lastLine.split()
	mx = float(data[1])
	my = float(data[2])
	mz = float(data[3])
	X = 0.600566
	Y = -0.00193843
	Z = -0.799573
	distance = sqrt((X-mx)**2+(Y-my)**2+(Z-mz)**2)
	threshold = 1e-5
	success = distance < threshold
	print("feeLLGood terminated successfully")
	print(f"final average reduced magnetization = ({mx}, {my}, {mz})")
	print("distance from expected value        = %.2e" % distance)
	print("threshold of acceptability          = %.2e" % threshold)
else:
	print("feeLLGood failed: exit status = " + str(val.returncode))
	success = False

# Use ANSI colors if printing to a terminal.
if sys.stdout.isatty():
	green  = "\x1b[1;32m"
	red    = "\x1b[1;31m"
	normal = "\x1b[m"
else:
	green  = ""
	red    = ""
	normal = ""

# Report success status.
if success:
	print(f"test {green}PASSED{normal}")
	sys.exit(0)
else:
	print(f"test {red}FAILED{normal}")
	sys.exit(1)
