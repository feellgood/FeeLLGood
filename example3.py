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
[0,phi,1], [0,phi,-1], [0,-phi,1], [0,-phi,-1]])

r1 = [ pts[0], pts[2], pts[3], pts[1] ] 
r2 = [ pts[4], pts[6], pts[7], pts[5] ]
r3 = [ pts[8], pts[10], pts[11], pts[9] ]

p3D = [r1,r2,r3]
#print('p3D=',p3D)
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


fc = [(1,0,0,0.2),(0,1,0,0.2),(0,0,1,0.2)]
ax.add_collection3d( Poly3DCollection( p3D, facecolors=fc, linewidths=1))
plt.show()

