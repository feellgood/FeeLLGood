import sys
import subprocess
import json
import numpy as np

sys.path.insert(0,'../tools')
from meshMaker import Cylinder
from settingsMaker import Settings


mesh_size = 0.5
r = 1
t=10


meshFileName = 'my_beautiful_cylinder.msh'
vol_regName = "bulk"
surf_regName = "frontier(bulk)"

cyl = Cylinder(r,t,mesh_size,surf_regName,vol_regName)

#surfaces for boundary conditions
bc1_regName = "S_left"
bc2_regName = "S_right"
cyl.addEdgeSurf(bc1_regName,bc2_regName)


cyl.make(meshFileName)

settings = {
    "outputs": {
        "file_basename": "test_ex5",
        "evol_time_step": 1e-12,
        "final_time": 2e-11,
        "evol_columns": [ "t", "<Mx>", "<My>", "<Mz>", "E_tot" ],
        "take_photo": False
    },
    "mesh": {
        "filename": meshFileName,
        "scaling_factor": 1e-9,
        "volume_regions": { vol_regName: { "alpha_LLG" : 0.02 } },
        "surface_regions": { surf_regName: {} }
    },
    "initial_magnetization": ["x", "y", 0],
    "Bext": [1, 0, -1],
    "finite_element_solver": { "nb_threads": 32 },
    "demagnetizing_field_solver": { "nb_threads": 32 },
    "time_integration": {
        "min(dt)": 5e-14,
        "max(dt)": 3.5e-13,
        "max(du)": 0.1
    },
    "spin_transfer_torque" : {
    "enable": True,
    "sigma": 1,
    "dens_state": 42,
    "beta": 1,
    "l_J": 1,
    "l_sf": 1,
    "V_file": True,
    "boundary_conditions": { bc1_regName: 1.2345, bc2_regName: -3.14 }
    }
}


#with open("rect_stt.json",'w') as outfile:
#    json.dump(settings,outfile,indent = 4)

val = subprocess.run(["../feellgood","-v", "-"], input=json.dumps(settings), text=True)

