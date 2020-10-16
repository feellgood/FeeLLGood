# python -m pip install -U pip
# sudo python -m pip install -U matplotlib

from math import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

# https://fr.wikipedia.org/wiki/Icosa%C3%A8dre
phi = (1+sqrt(5))/2
pts = np.array([[phi,1,0], [phi,-1,0], [-phi,1,0], [-phi,-1,0],
[1,0,phi], [1,0,-phi], [-1,0,phi], [-1,0,-phi],
[0,phi,1], [0,phi,-1], [0,-phi,1], [0,-phi,-1],
[0,0,0]] )

f_idx = [ [0,8,4], [0,4,1], [0,1,5], [0,5,9], [0,9,8], [8,9,2], [8,2,6], [8,6,4], [6,10,4], [6,3,10], [6,2,3], [10,3,11], [10,11,1], [10,1,4], [1,11,5], [11,3,7], [11,7,5], [3,2,7], [7,9,2],  [5,7,9] ]

r1 = [ pts[0], pts[2], pts[3], pts[1] ] 
r2 = [ pts[4], pts[5], pts[7], pts[6] ]
r3 = [ pts[8], pts[10], pts[11], pts[9] ]

triangle3D = []
for i in range(0,len(f_idx)):
	triangle3D.append( [pts[ f_idx[i][0] ],pts[ f_idx[i][1] ],pts[ f_idx[i][2] ]] )

meshFile = open("icosahedron.msh",'w')
meshFile.write("$MeshFormat\n2.2\t0\t8\n$EndMeshFormat\n")
meshFile.write("$Nodes\n")
meshFile.write(str(len(pts))+"\n")
for i in range(0,len(pts)):
	meshFile.write( str(i+1)+"\t"+str(pts[i][0])+"\t"+str(pts[i][1])+"\t"+str(pts[i][2]) + "\n" )
meshFile.write("$EndNodes\n")
meshFile.write("$Elements\n")
meshFile.write(str(2*len(f_idx))+"\n")
for i in range(0,len(f_idx)):
	meshFile.write( str(i+1)+"\t4\t2\t300\t1\t13\t"+str(1+f_idx[i][0])+"\t"+str(1+f_idx[i][1])+"\t"+str(1+f_idx[i][2]) + "\n" )
for i in range(0,len(f_idx)):
	meshFile.write( str(len(f_idx)+i+1)+"\t2\t2\t200\t1\t"+str(1+f_idx[i][0])+"\t"+str(1+f_idx[i][1])+"\t"+str(1+f_idx[i][2]) + "\n" )

meshFile.write("$EndElements\n")
meshFile.close()



fig = plt.figure()
#ax = plt.axes(projection='3d')
#x = []
#y = []
#z = []
#for i in range(0,len(pts)):
#	x.append(pts[i][0])
#	y.append(pts[i][1])
#	z.append(pts[i][2])
#ax.scatter(x,y,z,color='r')
#print(pts)
#print('x=',x)
#print('y=',y)
#print('z=',z)

ax = fig.add_subplot(111,projection='3d')

ax.set_xlabel('X')
ax.set_xlim3d(-2, 2)
ax.set_ylabel('Y')
ax.set_ylim3d(-2, 2)
ax.set_zlabel('Z')
ax.set_zlim3d(-2, 2)

fc = []
for i in range(1,21):
	fc.append( (i/20,0,0,0.2) )

triangle3D.append(r1)
triangle3D.append(r2)
triangle3D.append(r3)

fc.append( (1,1,0,0.4) )
fc.append( (0,1,0,0.4) )
fc.append( (0,0,1,0.4) )
#print(fc)
ax.add_collection3d( Poly3DCollection( triangle3D, facecolors=fc, linewidths=1))
plt.show()

