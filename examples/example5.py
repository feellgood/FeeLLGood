import sys
import subprocess
import json
import numpy as np

sys.path.insert(0,'../tools')
from cuboidMaker import Cuboid
from settingsMaker import Settings

rectangle = Cuboid([-2,-2,-2],[2,2,2],4,4,4)
meshFileName = 'rectangle.msh'
vol_regName = "300"
surf_regName = "200"

bc1_regName = "201"
bc2_regName = "202"

settings = {
    "outputs": {
        "file_basename": "test_ex5",
        "evol_time_step": 1e-12,
        "evol_columns": [ "t", "<Mx>", "<My>", "<Mz>", "E_tot" ],
        "take_photo": False
    },
    "mesh": {
        "filename": meshFileName,
        "scaling_factor": 1e-9,
        "volume_regions": { "300": { "alpha_LLG" : 0.02 } },
        "surface_regions": { "200": {} }
    },
    "initial_magnetization": ["x", "y", 0],
    "Bext": [1, 0, -1],
    "finite_element_solver": { "nb_threads": 32 },
    "demagnetizing_field_solver": { "nb_threads": 32 },
    "time_integration": {
        "final_time": 2e-11,
        "min(dt)": 5e-14,
        "max(dt)": 3.5e-13,
        "max(du)": 0.1
    },
    "spin_transfer_torque" : {
    "enable": True,
    "volume_region_reference": vol_regName,
    "sigma": 1,
    "dens_state": 42,
    "beta": 1,
    "l_J": 1,
    "l_sf": 1,
    "boundary_conditions": { bc1_regName: 1.2345, bc2_regName: -3.14 }
    }
}

rectangle.add_sub_surface(bc1_regName, lambda x,y,z: z== -2 )
rectangle.add_sub_surface(bc2_regName, lambda x,y,z: z== 2 )
rectangle.make(meshFileName,vol_regName,surf_regName)

#with open("rect_stt.json",'w') as outfile:
#    json.dump(settings,outfile,indent = 4)

val = subprocess.run(["../feellgood","-v", "-"], input=json.dumps(settings), text=True)

