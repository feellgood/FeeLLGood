#!/usr/bin/python3

import os
import sys
import json
import subprocess
from math import sqrt

# Move to the script's directory.
os.chdir(sys.path[0])

MaxNbThreads = int(subprocess.check_output(["getconf","_NPROCESSORS_ONLN"]))
nbThreads = (MaxNbThreads + 4) // 3 # ellipsoid.msh is small, maxnbthreads is counterproductive

settings = {
    "outputs": {
        "file_basename": "full_test",
        "evol_time_step": 1e-12,
        "final_time": 2e-11,
        "mag_config_every": False
    },
    "mesh": {
        "filename": "../examples/ellipsoid.msh",
        "length_unit": 1e-10,
        "volume_regions": { "ellipsoid_volume": {} },
        "surface_regions": { "ellipsoid_surface": {} }
    },
    "initial_magnetization": [0, 0, 1],
    "Bext": [1, 0, -1],
    "finite_element_solver": { "nb_threads": nbThreads },
    "demagnetizing_field_solver": { "nb_threads": nbThreads },
    "time_integration": {
        "min(dt)": 5e-14,
        "max(dt)": 3.5e-13,
        "max(du)": 0.1
    }
}

val = subprocess.run(["../feellgood", "--seed", "2", "-"], input=json.dumps(settings), text=True)

#(devNote) to avoid -fsanitize=leak ASLR bug, we turn off address randomizer
#val = subprocess.run(["setarch", "--addr-no-randomize", "../feellgood", "--seed", "2", "-"], input=json.dumps(settings), text=True)

if(val.returncode==0):
    with open("full_test.evol","r") as f:
        for line in f:
            pass
        lastLine = line
    f.close()
    data = lastLine.split()
    mx = float(data[1])
    my = float(data[2])
    mz = float(data[3])
    X = 0.323478
    Y = 0.377806
    Z = -0.867539
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
