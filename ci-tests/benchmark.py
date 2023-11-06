#!/usr/bin/python3

import os
import sys
import json
import subprocess
import timeit

from feellgood.meshMaker import Cylinder

def makeSettings(mesh, surface_name, volume_name, nbThreads):
    settings = {
        "outputs": {
            "file_basename": "benchmark",
            "evol_time_step": 1e-12,
            "final_time": 1e-10,
            "mag_config_every": False
        },
        "mesh": {
            "filename": mesh,
            "length_unit": 1e-9, # we use nanometers
            "volume_regions": { volume_name: {} },
            "surface_regions": { surface_name: {} }
        },
        "initial_magnetization": [0, 0, 1],
        "Bext": [1, 0, 1],
        "finite_element_solver": { "nb_threads": nbThreads },
        "demagnetizing_field_solver": { "nb_threads": nbThreads },
        "time_integration": {
            "min(dt)": 5e-18,
            "max(dt)": 1e-12,
            "max(du)": 0.1
        }
    }
    return settings

def task2test(settings):
    """ feellgood runs in a subprocess with seed=2 """
    val = subprocess.run(["../feellgood", "--seed", "2", "-"], input=json.dumps(settings), text=True)
    return val

def bench(outputFileName, elt_sizes, listNbThreads):
    meshFileName = "cylinder.msh"
    surface_name = "surface"
    volume_name = "volume"
    height = 32
    radius = 3.0 * height
    with open(outputFileName, 'w') as f:
        for elt_size in elt_sizes:
            f.write(str(elt_size) + '\t')
            mesh = Cylinder(radius, height, elt_size, surface_name, volume_name)
            mesh.make(meshFileName)
            for nbThreads in listNbThreads:
                settings = makeSettings(meshFileName, surface_name, volume_name, nbThreads)
                t = timeit.timeit(lambda: task2test(settings), number=1)
                if nbThreads == listNbThreads[-1]:
                    f.write(str(t) + '\n')
                else:
                    f.write(str(t) + '\t')
        f.close()

if __name__ == '__main__':
    os.chdir(sys.path[0])

    elt_sizes = [2.5, 3.0, 3.5, 4.0]
    listNbThreads = [1, 2, 4, 8, 16, 32, 64]
    bench('benchmark.txt', elt_sizes, listNbThreads)
