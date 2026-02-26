#!/usr/bin/python3

import os
import sys
import json
import subprocess
from math import sqrt

# Distance between two vectors of R^3.
def distance_R3(v1, v2):
    return sqrt((v2[0] - v1[0])**2 + (v2[1] - v1[1])**2 + (v2[2] - v1[2])**2)

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
    expected = { "m": [0.307432, 0.476202, -0.823843], "E_tot": -4.446980e-19 }
else:
    expected = { "m": [0.307542, 0.475877, -0.823989], "E_tot": -4.445605e-19 }

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
    data = list(map(float, lastLine.split()))
    m = [data[1], data[2], data[3]]
    E_tot = data[7]
    m_error = distance_R3(expected["m"], m)
    e_error = abs((E_tot - expected["E_tot"]) / expected["E_tot"])  # relative error
    threshold = 1e-5
    success = m_error < threshold and e_error < threshold
    print("feeLLGood terminated successfully")
    print(f"final average reduced magnetization = ({m[0]:.6f}, {m[1]:.6f}, {m[2]:.6f})")
    print(f"final total energy                  = {E_tot:.6e}")
    print(f"error on magnetization              = {m_error:.2e}")
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
