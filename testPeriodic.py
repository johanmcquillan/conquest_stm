
import numpy as np
from packages.vector import Vector
from packages.atomic import Atom, Radial
from packages.smart_dict import SmartDict
from copy import deepcopy

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
from mpl_toolkits.mplot3d import Axes3D

length = 20

v = Vector(length, length, length)

def constrain_vector(vector):
	"""Return a vector that is constrained within simulation cell"""
	new_vector = deepcopy(vector)

	# Check if vector components are greater than half of cell sides
	# If greater, add or subtract cell length

	if vector.x > v.x/2:
		new_vector -= v.project_x()
	elif vector.x <= -v.x/2:
		new_vector += v.project_x()

	if vector.y > v.y/2:
		new_vector -= v.project_y()
	elif vector.y <= -v.y/2:
		new_vector += v.project_y()

	if vector.z > v.z/2:
		new_vector -= v.project_z()
	elif vector.z <= -v.z/2:
		new_vector += v.project_z()

	# If vector is unchanged return
	# If vector has changed, constrain
	if new_vector == vector:
		return new_vector
	else:
		print 'Changing ', vector
		return constrain_vector(new_vector)

s = Vector(-21, -22, 53)
print constrain_vector(s)

xMesh, yMesh, zMesh = np.mgrid[0:length:1, 0:length:1, 0:length:1]
boolMesh = np.empty_like(xMesh, dtype=bool)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

radial = SmartDict()
radial[0][1] = Radial(1, 0, 1, [0, 1, 2, 3, 4, 5], [1, 5, 10, 3, 1, 0], 5)

print 'Making Atom'

atom = Atom('test', Vector(10, 10, 10), radial)

print 'Made Atom'

for i in range(length):
	for j in range(length):
		for k in range(length):
			x = xMesh[i, j, k]
			y = yMesh[i, j, k]
			z = zMesh[i, j, k]
			r = Vector(x, y, z)
			relative_r = constrain_vector(r - atom.atom_pos)
			if atom.within_cutoff_relative(relative_r):
				ax.scatter(x, y, z, 'o')

ax.set_xlim3d(0, length)
ax.set_ylim3d(0, length)
ax.set_zlim3d(0, length)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.show()
