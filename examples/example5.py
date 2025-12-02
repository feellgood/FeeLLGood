#!/usr/bin/env python3
import sys
import subprocess
import json
import numpy as np
import gmsh

verboseGmsh = False

justMagnet = False #if False magnet + electrode is meshed, else just the magnet, json also modified

#all dimensions are in nm
mesh_size = 15

class cylinder:
    def __init__(self,n,h,r):
        '''
        define a geometric cylinder of radius r and length h, named n (string) 
        '''
        self.name = n # region volume name
        self.height = h
        self.radius = r

def buildCircle(r,z0,orientation = True):
    '''
    create four arc circles composing a circle zero centered in the plane z=z0
    if orientation is True, return positive indices else negative (to be able to create holes with indirect orientation)
    '''
    p_origin = gmsh.model.geo.addPoint(0,0,z0,mesh_size)
    start_pt = gmsh.model.geo.addPoint(r,0,z0,mesh_size)
    interim_pt = gmsh.model.geo.addPoint(0,r,z0,mesh_size)
    interim1_pt = gmsh.model.geo.addPoint(-r,0,z0,mesh_size)
    interim2_pt = gmsh.model.geo.addPoint(0,-r,z0,mesh_size)
    arc0 = gmsh.model.geo.addCircleArc(start_pt,p_origin,interim_pt)
    arc1 = gmsh.model.geo.addCircleArc(interim_pt,p_origin,interim1_pt)
    arc2 = gmsh.model.geo.addCircleArc(interim1_pt,p_origin,interim2_pt)
    arc3 = gmsh.model.geo.addCircleArc(interim2_pt,p_origin,start_pt)
    gmsh.model.geo.remove([(p_origin,0)])
    if orientation:
        return([arc0,arc1,arc2,arc3])
    else:
        return([-arc3,-arc2,-arc1,-arc0])

def buildDisk(circle):
    '''
    create a disk from a circle
    '''
    curvedLoop = gmsh.model.geo.addCurveLoop(circle) # curvedLoop is an index (integer)
    surf = gmsh.model.geo.addPlaneSurface([curvedLoop]) # surf is an index (integer)
    return surf

meshFileName = "nanowire.msh"
# dimensions of the whole nanopillar = nanowire + electrode
l_n = 100
r_n = 20
nw = cylinder("magnet",l_n,r_n) #nanowire
l_e = 40
r_e = 30
e = cylinder("metal",l_e,r_e) #electrode
if r_n > r_e:
    print("Error: diameter of the magnetic nanopillar exceed electrode diameter.")
    exit(1)

# cross section of the mesh (xOz) plane:
#
#      ---  z= nw.height (#202) boundary condition J = 1.0
#     |   |
#     |   |
#     |   |    magnet (z>0) radius = r_n
#     |   |
#     |---|  z=0
#     |   |    electrode (z<0) radius = r_e
#      ---   z= -e.height (#212)  boundary condition V = 0.0

gmsh.initialize()
gmsh.option.setNumber("General.Terminal",verboseGmsh)
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.model.add("nanowireWithElectrode")

z0=0
circle0 = buildCircle(nw.radius,z0)
surf = buildDisk(circle0)
out = gmsh.model.geo.extrude([(2,surf)],0,0,nw.height) # 2 is the dimension of the object refered by index surf
gmsh.model.geo.synchronize() # we have to sync before calling addPhysicalGroup

surface_tag = 200 # this surface is frontier(magnetic volume)
surf_regName = "frontier(magnet)"
gmsh.model.addPhysicalGroup(2,[surf, out[0][1], out[2][1], out[3][1], out[4][1], out[5][1]],surface_tag)
gmsh.model.setPhysicalName(2,surface_tag,surf_regName)

if not justMagnet:
    surface_tag = 202 # should be disk in z=nw.height plane : boundary condition for electrostatic problem
    surface_top_name = "top_mag_pillar"
    gmsh.model.addPhysicalGroup(2,[ out[0][1] ],surface_tag)
    gmsh.model.setPhysicalName(2,surface_tag,surface_top_name)

volume_tag = 300 # magnetic volume
gmsh.model.addPhysicalGroup(3,[out[1][1]],volume_tag)
gmsh.model.setPhysicalName(3,volume_tag,nw.name)

if r_n < r_e:
    circle1 = list(map(lambda x: -x,circle0))  # change all sign of circle0 to define the hole
    circle2 = buildCircle(e.radius,z0)
    curvedLoop = gmsh.model.geo.addCurveLoop(circle2)
    hole = gmsh.model.geo.addCurveLoop(circle1)
    surf2 = gmsh.model.geo.addPlaneSurface([curvedLoop,hole])
    out2 = gmsh.model.geo.extrude([(1,circle2[0]),(1,circle2[1]),(1,circle2[2]),(1,circle2[3])],0,0,-e.height)
    gmsh.model.geo.synchronize()
    circle4 = buildCircle(e.radius,z0-e.height)
    surf3 = buildDisk(circle4)
    gmsh.model.geo.synchronize() # we have to sync before calling addPhysicalGroup

    frontier_metal = [surf,surf2,out2[1][1],out2[5][1],out2[9][1],out2[13][1],surf3]
    surf_loop = gmsh.model.geo.addSurfaceLoop(frontier_metal)
    gmsh.model.geo.removeAllDuplicates() # we have to, otherwise some error occurs while 3D meshing, because of "intersection of line and point", weird...
    metal_vol = gmsh.model.geo.addVolume([surf_loop])
    gmsh.model.geo.synchronize()
elif r_n == r_e:
    out2 = gmsh.model.geo.extrude([(2,surf)],0,0,-e.height) # 2 is the dimension of the object refered by index surf
    gmsh.model.geo.synchronize() # we have to sync before calling addPhysicalGroup
    metal_vol = out2[1][1]
    surf3 = out2[0][1]

if not justMagnet:
    surface_tag = 212 # disk in z= z0-e.height plane : boundary condition for electrostatic problem
    surface_top_name2 = "metal_electrode"
    gmsh.model.addPhysicalGroup(2,[surf3],surface_tag)
    gmsh.model.setPhysicalName(2,surface_tag,surface_top_name2)

    volume_tag = 310 # normal metal volume
    gmsh.model.addPhysicalGroup(3,[metal_vol],volume_tag)
    gmsh.model.setPhysicalName(3,volume_tag,e.name)

gmsh.model.geo.synchronize() # we have to synchronize before the call to 'generate' to build the mesh
gmsh.model.mesh.generate(3) # 3 is the dimension of the mesh
        
gmsh.write(meshFileName)
print(f"Generated {meshFileName}")
gmsh.finalize()

# warning: Cobalt l_sd was randomly choosen, just to have spin diffusion solver converges
# polarization vector P is modelling a polarizer "outside" the mesh, such as a polarizing layer. It
# is attributed to the adjacency surface of the polarizing layer and the magnet.
settings = {
    "outputs": {
        "file_basename": "test_ex5",
        "evol_time_step": 1e-13,
        "final_time": 1e-11,
        "evol_columns": [ "t", "<Mx>", "<My>", "<Mz>", "E_tot" ],
        "mag_config_every": 20
    },
    "mesh": {
        "filename": meshFileName,
        "length_unit": 1e-9,
        "volume_regions": {
            nw.name: { "Ae": 1e-11, "Js": 1, "alpha_LLG": 0.05, "P": 0.7, "sigma": 1.7e7,
                       "dens_state": 1.34e47, "l_sd": 10e-9, "l_sf": 12.5e-9} #Co
            },
        "surface_regions": {
            surf_regName: {}
            }
    },
    "initial_magnetization": ["-y", "x", "2e-9"],
    "Bext": [0, 0, 0],
    "demagnetizing_field_solver": { "nb_threads": 32 },
    "time_integration": {
        "min(dt)": 0.5e-16,
        "max(dt)": 1e-11
    },
    "spin_accumulation": { "enable": True, "V_file": True }
}

if not justMagnet:
    settings["mesh"]["volume_regions"][e.name] = { "Ae": 0, "Js": 0, "P": 0.7, "sigma":
    5.8e7, "dens_state": 1.34e74, "l_sf": 350e-9 } #Cu
    settings["mesh"]["surface_regions"][surface_top_name] = { "J": 1.0e12, "P":[0,1,0] }
    settings["mesh"]["surface_regions"][surface_top_name2] = { "V": 0.0, "s": [0,0,0] }

jsonFileName = "nanowire_spinAcc.json"
with open(jsonFileName,'w') as outfile:
    json.dump(settings,outfile,indent = 4)

print(f"Generated {jsonFileName}")

