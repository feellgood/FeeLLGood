#!/usr/bin/env python3
import json
from feellgood.meshMaker import Tube

r1= 50 #inner radius
r2= 70 #outer radius
l= 200 #length
meshSize= 4
meshFileName = "tube.msh"
volRegionName = "bulk tube"
surfRegionName = "frontier(" + volRegionName + ")"
t = Tube(r1,r2,l,meshSize,surfRegionName,volRegionName)
#t.addEdgeSurf("leftS","rightS")
t.make(meshFileName)
settings = {
    "outputs": {
        "file_basename": "tube",
        "evol_time_step": 1e-13,
        "final_time": 1e-11,
        "evol_columns": [ "t", "<Mx>", "<My>", "<Mz>", "E_tot" ],
        "mag_config_every": 20
    },
    "mesh": {
        "filename": meshFileName,
        "length_unit": 1e-9,
        "volume_regions": {
            volRegionName: { "Ae": 1e-11, "Ms": 800e3, "alpha_LLG": 0.05 }
            },
        "surface_regions": {
            surfRegionName: {}
            }
    },
    "initial_magnetization": ["-y*sign(z)", "x*sign(z)", "2e-9"],
    "Bext": [0, 0, 0],
    "demagnetizing_field_solver": { "nb_threads": 32 },
    "time_integration": {
        "min(dt)": 0.5e-16,
        "max(dt)": 1e-11
    }
}
jsonFileName = "tube.json"
with open(jsonFileName,'w') as outfile:
    json.dump(settings,outfile,indent = 4)
print(f"Generated {jsonFileName}")

