#!/usr/bin/env python3
import sys
import subprocess
import json
import numpy as np

from feellgood.meshMaker import Cylinder

#nanowire dimensions in nm
mesh_size = 4
height = 100
radius = 20

meshFileName = "nanowire.msh"
vol_regName = "bulk"
surf_regName = "frontier(bulk)"

cyl = Cylinder(radius, height, mesh_size, surf_regName, vol_regName)

#surfaces for boundary conditions
surface_base_name = "ground"
surface_top_name = "top_electrode"

#create cylinder edge surfaces as 2D geometric entities in the mesh, to define boundary conditions
cyl.addEdgeSurf(surface_base_name,surface_top_name)

cyl.make(meshFileName)

settings = {
    "outputs": {
        "file_basename": "test_ex5",
        "evol_time_step": 1e-12,
        "final_time": 2e-11,
        "evol_columns": [ "t", "<Mx>", "<My>", "<Mz>", "E_tot" ],
        "mag_config_every": False
    },
    "mesh": {
        "filename": meshFileName,
        "length_unit": 1e-9,
        "volume_regions": { vol_regName: { "Ae": 1e-11, "Js": 1, "alpha_LLG" : 0.05 } },
        "surface_regions": { surf_regName: {}, surface_base_name:{"V":0.0}, surface_top_name:{"J":1.1} }
    },
    "initial_magnetization": [1, 0, 1],
    "Bext": [0, 0, 0],
    "finite_element_solver": { "nb_threads": 32 },
    "demagnetizing_field_solver": { "nb_threads": 32 },
    "time_integration": {
        "min(dt)": 0.5e-16,
        "max(dt)": 1e-11
    },
    "spin_transfer_torque" : {
        "enable": True,
        "sigma": 1,
        "dens_state": 42,
        "beta": 1,
        "l_J": 1,
        "l_sf": 1,
        "V_file": False
    }
}

#with open("rect_stt.json",'w') as outfile:
#    json.dump(settings,outfile,indent = 4)

val = subprocess.run(["../feellgood","-v", "-"], input=json.dumps(settings), text=True)

