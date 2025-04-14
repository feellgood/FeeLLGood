from math import sqrt
import numpy as np

from feellgood.settingsMaker import Settings

phi = (1+sqrt(5))/2
pts = np.array([[phi,1,0], [phi,-1,0], [-phi,1,0], [-phi,-1,0],[1,0,phi], [1,0,-phi], [-1,0,phi], [-1,0,-phi],[0,phi,1], [0,phi,-1], [0,-phi,1], [0,-phi,-1], [0,0,0]] )

f_idx = [ [0,8,4], [0,4,1], [0,1,5], [0,5,9], [0,9,8], [8,9,2], [8,2,6], [8,6,4], [6,10,4], [6,3,10], [6,2,3], [10,3,11], [10,11,1], [10,1,4], [1,11,5], [11,3,7], [11,7,5], [3,2,7], [7,9,2],  [5,7,9] ]

pts2 = pts + np.array([4.0,0,0])

meshFileName = "twinIco.msh"
volRegionName = "bulk"
surfRegionName = "frontier(bulk)"
volRegionName2 = "bulk2"
surfRegionName2 = "frontier(bulk2)"

meshFile = open(meshFileName,'w')
meshFile.write("$MeshFormat\n2.2\t0\t8\n$EndMeshFormat\n")
meshFile.write("$PhysicalNames\n4\n")
surfRegionTag = 200
surfRegionTag2 = 201
volRegionTag = 300
volRegionTag2 = 301
meshFile.write("2 " + str(surfRegionTag) + ' \"' + surfRegionName + '\"\n')
meshFile.write("2 " + str(surfRegionTag2) + ' \"' + surfRegionName2 + '\"\n')
meshFile.write("3 " + str(volRegionTag) + ' \"' + volRegionName + '\"\n')
meshFile.write("3 " + str(volRegionTag2) + ' \"' + volRegionName2 + '\"\n')
meshFile.write("$EndPhysicalNames\n")
meshFile.write("$Nodes\n" + str(len(pts) + len(pts2))+"\n")
for i in range(0,len(pts)):
    meshFile.write( str(i+1)+"\t"+str(pts[i][0])+"\t"+str(pts[i][1])+"\t"+str(pts[i][2]) + "\n" )
for i in range(0,len(pts2)):
    meshFile.write( str(i+14)+"\t"+str(pts2[i][0])+"\t"+str(pts2[i][1])+"\t"+str(pts2[i][2]) + "\n" )
meshFile.write("$EndNodes\n$Elements\n")
meshFile.write(str(4*len(f_idx))+"\n")
for i in range(0,len(f_idx)):
     meshFile.write( str(i+1)+"\t4\t2\t" + str(volRegionTag) + "\t1\t13\t"+str(1+f_idx[i][0])+"\t"+str(1+f_idx[i][1])+"\t"+str(1+f_idx[i][2]) + "\n" )
offset = len(f_idx)
for i in range(0,len(f_idx)):
     meshFile.write( str(offset+i+1)+"\t4\t2\t" + str(volRegionTag2) + "\t2\t26\t"+str(14+f_idx[i][0])+"\t"+str(14+f_idx[i][1])+"\t"+str(14+f_idx[i][2]) + "\n" )
offset += len(f_idx)
for i in range(0,len(f_idx)):
     meshFile.write( str(offset+i+1)+"\t2\t2\t" + str(surfRegionTag) + "\t1\t"+str(1+f_idx[i][0])+"\t"+str(1+f_idx[i][1])+"\t"+str(1+f_idx[i][2]) + "\n" )
offset += len(f_idx)
for i in range(0,len(f_idx)):
     meshFile.write( str(offset+i+1)+"\t2\t2\t" + str(surfRegionTag2) + "\t2\t"+str(14+f_idx[i][0])+"\t"+str(14+f_idx[i][1])+"\t"+str(14+f_idx[i][2]) + "\n" )

meshFile.write("$EndElements\n")
meshFile.close()
print(f"Generated {meshFileName}: two icosaedrons")

mySettings = Settings(meshFileName)

mySettings["outputs"]["file_basename"] = "twinIco"

mySettings["outputs"]["evol_time_step"] = 1e-11
mySettings["outputs"]["final_time"] = 5e-9
mySettings["outputs"]["evol_columns"] = [ "t", "<Mx>", "<My>", "<Mz>", "E_tot" ]

mySettings["outputs"]["mag_config_every"] = False

mySettings.createVolRegion( volRegionName )
mySettings.createSurfRegion( surfRegionName )
mySettings["mesh"]["volume_regions"][volRegionName]["Js"] = 1.0
mySettings["mesh"]["volume_regions"][volRegionName]["Ae"] = 10e-12
mySettings["mesh"]["volume_regions"][volRegionName]["alpha_LLG"] = 0.1

mySettings.createVolRegion( volRegionName2 )
mySettings.createSurfRegion( surfRegionName2 )
mySettings["mesh"]["volume_regions"][volRegionName2]["Js"] = 0.5
mySettings["mesh"]["volume_regions"][volRegionName2]["Ae"] = 10e-12
mySettings["mesh"]["volume_regions"][volRegionName2]["alpha_LLG"] = 0.1

mySettings["time_integration"]["min(dt)"] = 5e-18
mySettings["time_integration"]["max(dt)"] = 5e-10

mySettings["initial_magnetization"] = [1, 0, 1]

mySettings.write('twinIco.json')

