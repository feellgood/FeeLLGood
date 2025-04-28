#!/usr/bin/env python3
import sys
import subprocess
import json
import numpy as np
import gmsh

verboseGmsh = False

#all dimensions are in nm
mesh_size = 5

class cylinder:
    def __init__(self,h,r):
        self.height = h
        self.radius = r
#nanowire
nw = cylinder(100,20)

#electrode
e = cylinder(40,20)

# cross section of the mesh (xOz) plane:
#
#      ---  z= nw.height (#202) boundary condition V = 3.0
#     |   |
#     |   |
#     |   |    magnet (z>0)
#     |   |
#     |---|  z=0
#     |   |    electrode (z<0)
#      ---   z= -e.height (#212)  boundary condition J = 0.1

meshFileName = "nanowire.msh"
vol_regName = "magnet"
surf_regName = "frontier(magnet)"

vol_regName2 = "metal"
surf_regName2 = "frontier(metal)"

surface_top_name = "top_ground"

surface_top_name2 = "metal_top"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal",verboseGmsh)
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.model.add("nanowireWithElectrode")

p_origin = gmsh.model.geo.addPoint(0,0,0,mesh_size)

start_pt = gmsh.model.geo.addPoint(nw.radius,0,0,mesh_size)
interim_pt = gmsh.model.geo.addPoint(0,nw.radius,0,mesh_size)
interim1_pt = gmsh.model.geo.addPoint(-nw.radius,0,0,mesh_size)
interim2_pt = gmsh.model.geo.addPoint(0,-nw.radius,0,mesh_size)
        
big_circle = gmsh.model.geo.addCircleArc(start_pt,p_origin,interim_pt)
big_circle1 = gmsh.model.geo.addCircleArc(interim_pt,p_origin,interim1_pt)
big_circle2 = gmsh.model.geo.addCircleArc(interim1_pt,p_origin,interim2_pt)
big_circle3 = gmsh.model.geo.addCircleArc(interim2_pt,p_origin,start_pt)

gmsh.model.geo.remove([(p_origin,0)])

curvedLoop = gmsh.model.geo.addCurveLoop([big_circle,big_circle1,big_circle2,big_circle3]) # curvedLoop is an index (integer)
surf = gmsh.model.geo.addPlaneSurface([curvedLoop]) # surf is an index (integer)
out = gmsh.model.geo.extrude([(2,surf)],0,0,nw.height) # 2 is the dimension of the object refered by index surf
gmsh.model.geo.synchronize() # we have to sync before calling addPhysicalGroup

surface_tag = 200 # this surface is frontier(magnetic volume)
gmsh.model.addPhysicalGroup(2,[surf, out[0][1], out[2][1], out[3][1], out[4][1], out[5][1]],surface_tag)
gmsh.model.setPhysicalName(2,surface_tag,surf_regName)

#create cylinder edge surfaces as 2D geometric entities in the mesh, to define boundary conditions
#surface_tag_left = 201
#gmsh.model.addPhysicalGroup(2,[surf],surface_tag_left)
#gmsh.model.setPhysicalName(2,surface_tag_left,surface_base_name)

surface_tag = 202 # should be disk in z=nw.height plane
gmsh.model.addPhysicalGroup(2,[ out[0][1] ],surface_tag)
gmsh.model.setPhysicalName(2,surface_tag,surface_top_name)

volume_tag = 300 # magnetic volume
gmsh.model.addPhysicalGroup(3,[out[1][1]],volume_tag)
gmsh.model.setPhysicalName(3,volume_tag,vol_regName)

out2 = gmsh.model.geo.extrude([(2,surf)],0,0,-e.height) # 2 is the dimension of the object refered by index surf

gmsh.model.geo.synchronize() # we have to sync before calling addPhysicalGroup

surface_tag = 212 # should be disk in z= -e.height plane
gmsh.model.addPhysicalGroup(2,[ out2[0][1] ],surface_tag)
gmsh.model.setPhysicalName(2,surface_tag,surface_top_name2)

volume_tag = 310 # normal metal volume
gmsh.model.addPhysicalGroup(3,[out2[1][1]],volume_tag)
gmsh.model.setPhysicalName(3,volume_tag,vol_regName2)

gmsh.model.geo.synchronize() # we have to synchronize before the call to 'generate' to build the mesh
gmsh.model.mesh.generate(3) # 3 is the dimension of the mesh
        
gmsh.write(meshFileName)

# uncomment next line to see a graphic rendering of the mesh
#gmsh.fltk.run()

print(f"Generated {meshFileName}")
gmsh.finalize()

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
        "volume_regions": {
            vol_regName: { "Ae": 1e-11, "Js": 1, "alpha_LLG": 0.05, "beta": 0.7, "sigma": 1.7e7 }, #Co
            vol_regName2: { "Ae": 0, "Js": 0, "beta": 0.7, "sigma": 5.8e7 } #Cu
            },
        "surface_regions": {
            surf_regName: {},
            surface_top_name2:{ "J": 0.1 },
            surface_top_name:{ "V": 3.0 }
            }
    },
    "initial_magnetization": [1, 0, 1],
    "Bext": [0, 0, 0],
    "demagnetizing_field_solver": { "nb_threads": 32 },
    "time_integration": {
        "min(dt)": 0.5e-16,
        "max(dt)": 1e-11
    },
    "spin_accumulation": { "enable": True, "V_file": True }
}

jsonFileName = "nanowire_spinAcc.json"
with open(jsonFileName,'w') as outfile:
    json.dump(settings,outfile,indent = 4)

print(f"Generated {jsonFileName}")
#val = subprocess.run(["../feellgood","-v", "-"], input=json.dumps(settings), text=True)

