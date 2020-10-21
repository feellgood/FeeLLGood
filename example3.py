# python -m pip install -U pip
# sudo python -m pip install -U matplotlib

from math import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from settingsMaker import Settings

meshFileName = 'icosahedron.msh'

mySettings = Settings(meshFileName)

mySettings["outputs"]["file_basename"] = "ico"
mySettings["outputs"]["evol_time_step"] = 0.2e-8 
mySettings["outputs"]["evol_columns"] = [ "iter", "t", "dt", "<Mx>", "<My>", "<Mz>", "E_ex", "E_aniso", "E_demag", "E_tot" ]

#carefull here, volume and region names must be stringified integers
vol_region_name = "300"
surf_region_name = "200"

mySettings.createVolRegion( vol_region_name )
mySettings.createSurfRegion( surf_region_name )

# cubic anisotropy for volume region "300" (Fe : K3 = 4.2e4 J/m^3)
mySettings["mesh"]["volume_regions"]["300"]["K3"] = 4.2e4

#cosine directions [alpha,beta,gamma]
mySettings["mesh"]["volume_regions"]["300"]["alpha"] = [1.0,0.0,0.0]
mySettings["mesh"]["volume_regions"]["300"]["beta"] = [0.0,1.0,0.0]
mySettings["mesh"]["volume_regions"]["300"]["gamma"] = [0.0,0.0,1.0]

mySettings["initial_magnetization"] = {"Mx":"1","My":"1","Mz":"1"}

mySettings.write('ico_test.json')

# https://fr.wikipedia.org/wiki/Icosa%C3%A8dre
phi = (1+sqrt(5))/2
pts = np.array([[phi,1,0], [phi,-1,0], [-phi,1,0], [-phi,-1,0],
[1,0,phi], [1,0,-phi], [-1,0,phi], [-1,0,-phi],
[0,phi,1], [0,phi,-1], [0,-phi,1], [0,-phi,-1],
[0,0,0]] )

f_idx = [ [0,8,4], [0,4,1], [0,1,5], [0,5,9], [0,9,8], [8,9,2], [8,2,6], [8,6,4], [6,10,4], [6,3,10], [6,2,3], [10,3,11], [10,11,1], [10,1,4], [1,11,5], [11,3,7], [11,7,5], [3,2,7], [7,9,2],  [5,7,9] ]

meshFile = open(meshFileName,'w')
meshFile.write("$MeshFormat\n2.2\t0\t8\n$EndMeshFormat\n$Nodes\n")
meshFile.write(str(len(pts))+"\n")
for i in range(0,len(pts)):
	meshFile.write( str(i+1)+"\t"+str(pts[i][0])+"\t"+str(pts[i][1])+"\t"+str(pts[i][2]) + "\n" )
meshFile.write("$EndNodes\n$Elements\n")
meshFile.write(str(2*len(f_idx))+"\n")
for i in range(0,len(f_idx)):
	meshFile.write( str(i+1)+"\t4\t2\t" + vol_region_name + "\t1\t13\t"+str(1+f_idx[i][0])+"\t"+str(1+f_idx[i][1])+"\t"+str(1+f_idx[i][2]) + "\n" )
for i in range(0,len(f_idx)):
	meshFile.write( str(len(f_idx)+i+1)+"\t2\t2\t" + surf_region_name + "\t1\t"+str(1+f_idx[i][0])+"\t"+str(1+f_idx[i][1])+"\t"+str(1+f_idx[i][2]) + "\n" )

meshFile.write("$EndElements\n")
meshFile.close()

# the following part of the script is only to have a 3D view of the icosahedron, useless if you do not want to see it

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

ax.set_xlabel('X')
ax.set_xlim3d(-2, 2)
ax.set_ylabel('Y')
ax.set_ylim3d(-2, 2)
ax.set_zlabel('Z')
ax.set_zlim3d(-2, 2)

r1 = [ pts[0], pts[2], pts[3], pts[1] ] 
r2 = [ pts[4], pts[5], pts[7], pts[6] ]
r3 = [ pts[8], pts[10], pts[11], pts[9] ]

triangle3D = []
for i in range(0,len(f_idx)):
	triangle3D.append( [pts[ f_idx[i][0] ],pts[ f_idx[i][1] ],pts[ f_idx[i][2] ]] )

triangle3D.append(r1)
triangle3D.append(r2)
triangle3D.append(r3)

fc = []
for i in range(1,len(f_idx)+1):
	fc.append( (i/20,0,0,0.2) )

fc.append( (1,1,0,0.4) )
fc.append( (0,1,0,0.4) )
fc.append( (0,0,1,0.4) )

ax.add_collection3d( Poly3DCollection( triangle3D, facecolors=fc, linewidths=1))
plt.show()

