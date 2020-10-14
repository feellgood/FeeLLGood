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

zero = np.array( [0,0,0] )

r1 = [ pts[0], pts[2], pts[3], pts[1] ] 
r2 = [ pts[4], pts[5], pts[7], pts[6] ]
r3 = [ pts[8], pts[10], pts[11], pts[9] ]

triangle3D = []
triangle3D.append( [pts[0], pts[8], pts[4] ] )
triangle3D.append( [pts[0], pts[4], pts[1] ] )
triangle3D.append( [pts[0], pts[1], pts[5] ] )
triangle3D.append( [pts[0], pts[5], pts[9] ] )
triangle3D.append( [pts[0], pts[9], pts[8] ] )
triangle3D.append( [pts[8], pts[9], pts[2] ] )
triangle3D.append( [pts[8], pts[2], pts[6] ] )
triangle3D.append( [pts[8], pts[6], pts[4] ] )
triangle3D.append( [pts[6], pts[10], pts[4] ] )
triangle3D.append( [pts[6], pts[3], pts[10] ] )
triangle3D.append( [pts[6], pts[2], pts[3] ] )
triangle3D.append( [pts[10], pts[3], pts[11]] )
triangle3D.append( [pts[10], pts[11], pts[1]] )
triangle3D.append( [pts[10], pts[1], pts[4]] )
triangle3D.append( [pts[1], pts[11], pts[5]] )
triangle3D.append( [pts[11], pts[3], pts[7]] )
triangle3D.append( [pts[11], pts[7], pts[5]] )
triangle3D.append( [pts[3], pts[2], pts[7]] )
triangle3D.append( [pts[3], pts[6], pts[2]] )
triangle3D.append( [pts[5], pts[7], pts[9]] )

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

#fc.append( (1,1,0,0.4) )
#fc.append( (0,1,0,0.4) )
#fc.append( (0,0,1,0.4) )
#print(fc)
ax.add_collection3d( Poly3DCollection( triangle3D, facecolors=fc, linewidths=1))
plt.show()

