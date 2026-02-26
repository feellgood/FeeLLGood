#!/usr/bin/python3

import os
import sys
import json
import subprocess
from math import sqrt

# Move to the script's directory.
os.chdir(sys.path[0])
sys.path.insert(0, '../python-modules')

# Simulation settings. The number of threads depends on the number of available CPUs.
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
        "length_unit": 1e-9,
        "volume_regions": {
            "ellipsoid_volume": { "K": 3e5, "uk": [0, 1, 0] }
        },
        "surface_regions": { "ellipsoid_surface": {} }
    },
    "initial_magnetization": [0, 0, 1],
    "Bext": [1, 0, -1],
    "demagnetizing_field_solver": { "nb_threads": nbThreads },
    "time_integration": {
        "min(dt)": 5e-14,
        "max(dt)": 1e-12,
        "max(du)": 0.1
    }
}

# Expected outcome. It depends on whether feeLLGood was compiled in ONE_GAUSS_POINT mode.
val = subprocess.check_output(["../feellgood","--version"])
if b'ONE_GAUSS_POINT=ON' in val:
    X = 0.307432
    Y = 0.476202
    Z = -0.823843
    E = -4.446980e-19
else:
    X = 0.307542
    Y = 0.475877
    Z = -0.823989
    E = -4.445605e-19

# Run the simulation.
sys.stdout.flush()
val = subprocess.run(["../feellgood", "--seed", "2", "-"], input=json.dumps(settings), text=True)

#(devNote) to avoid -fsanitize=leak ASLR bug, we turn off address randomizer
#val = subprocess.run(["setarch", "--addr-no-randomize", "../feellgood", "--seed", "2", "-"], input=json.dumps(settings), text=True)

# Report results.
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
    e_tot = float(data[7])
    distance = sqrt((X-mx)**2+(Y-my)**2+(Z-mz)**2)
    e_error = abs((e_tot - E) / E)  # relative error
    threshold = 1e-5
    success = distance < threshold and e_error < threshold
    print("feeLLGood terminated successfully")
    print(f"final average reduced magnetization = ({mx:.6f}, {my:.6f}, {mz:.6f})")
    print(f"distance from expected value        = {distance:.2e}")
    print(f"final total energy                  = {e_tot:.6e}")
    print(f"relative error on total energy      = {e_error:.2e}")
    print(f"threshold of acceptability          = {threshold:.2e}")
else:
    print("feeLLGood failed: exit status = " + str(val.returncode))
    success = False

# Use ANSI colors if printing to a terminal.
if sys.stdout.isatty() and not os.getenv("NO_COLOR"):
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
