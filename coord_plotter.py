
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = 'x'
y = 'y'
z = 'z'
axes = [x, y, z]

latticePoints = {x : [.0, .0, .5, .5], y : [.0, .5, .0, .5], z : [.0, .5, .5, .0]}
basisPoints = {x : [], y : [], z : []}
basisVector = {x : .25, y: .25, z : .25}

latticePoints = {x : [.0, .0, .5, .5], y : [.0, .5, .0, .5], z : [.0, .5, .5, .0]}
basisPoints = {x : [], y : [], z : []}
basisVector = {x : .25, y: .25, z : .25}

for e in axes:
	for i in range(len(latticePoints[e])):
		basisPoints[e].append(latticePoints[e][i] + basisVector[e])

plot1 = False

if plot1:
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(latticePoints[x], latticePoints[y], latticePoints[z], c='r', marker='o')
	ax.scatter(basisPoints[x], basisPoints[y], basisPoints[z], c='b', marker='^')

else:


	layers = 2
	squares = 2
	gap = 0
	zFactor = 1.0 / (layers + 2*gap)
	xyFactor = 1.0 / squares
	surfacePoints = {x : [], y : [], z : []}

	for e in axes:
		for i in range(0, layers):
			for j in range(0, squares):
				for k in range(0, squares):
					for h in range(len(latticePoints[e])):
						if e == 'z':
							surfacePoints[e].append((gap + i + latticePoints[e][h])*zFactor)
							surfacePoints[e].append((gap + i + latticePoints[e][h] + basisVector[e])*zFactor)
						elif e == 'x':
							surfacePoints[e].append((j + latticePoints[e][h])*xyFactor)
							surfacePoints[e].append((j + latticePoints[e][h] + basisVector[e])*xyFactor)
						elif e == 'y':
							surfacePoints[e].append((k + latticePoints[e][h])*xyFactor)
							surfacePoints[e].append((k + latticePoints[e][h] + basisVector[e])*xyFactor)



	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	bonds = True
	if bonds:
		for i in range(0, len(surfacePoints[x])):
			for j in range(0, i):
				x1 = surfacePoints[x][i]
				y1 = surfacePoints[y][i]
				z1 = surfacePoints[z][i]
				x2 = surfacePoints[x][j]
				y2 = surfacePoints[y][j]
				z2 = surfacePoints[z][j]

				xdistance = x2 - x1
				ydistance = y2 - y1
				zdistance = z2 - z1

				distance = np.sqrt(xdistance**2 + ydistance**2 + zdistance**2)
				if (i != j) and (distance < 0.0):
					if z1 > 0.8 or z2 > 0.8:
						ax.plot([x1, x2], [y1, y2], [z1, z2], color='b')
					else:
						ax.plot([x1, x2], [y1, y2], [z1, z2], color='g')
	cell = True
	if cell:
		lxy = [.25, .75]
		lz = [.0, .5]
		for p in [[.25, .25], [.75, .75]]:
			for q in [[.0, .0], [.5, .5]]:
				for e in ['x', 'y']:	
					if e == 'x':
						ax.plot(lxy, p, q, color='g')
					elif e == 'y':
						ax.plot(p, lxy, q, color='g')
				for p2 in [[.25, .25], [.75, .75]]:
					ax.plot(p, p2, lz, color='g')

		lxy = [.25, .75]
		lz = [.0, .5]
		for p in [[.25, .25], [.75, .75]]:
			for q in [[.0, .0], [.5, .5]]:
				for e in ['x', 'y']:	
					if e == 'x':
						ax.plot(lxy, p, q, color='g')
					elif e == 'y':
						ax.plot(p, lxy, q, color='g')
				for p2 in [[.25, .25], [.75, .75]]:
					ax.plot(p, p2, lz, color='g')

	for i in range(len(surfacePoints[x])):
		coordString = "( "+str(surfacePoints[x][i])+", "+str(str(surfacePoints[y][i]))+", "+str(str(surfacePoints[z][i]))+" )"
		print coordString

	ax.scatter(surfacePoints[x], surfacePoints[y], surfacePoints[z], c='r', marker='o')

ax.set_xlim3d(0, 1)
ax.set_ylim3d(0, 1)
ax.set_zlim3d(0, 1)


plt.show()
plt.close()
