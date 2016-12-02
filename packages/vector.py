import numpy as np


class Vector(object):
	"""3D cartesian vector class.

	Attributes:
		x (float): Component value in x direction
		y (float): Component value in y direction
		z (float): Component value in z direction
		r (float): Vector magnitude. Initially null - only calculates if needed, saves computational cost
	"""

	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
		self.r = None

	def subtract(self, V):
		dx = self.x - V.x
		dy = self.y - V.y
		dz = self.z - V.z
		return Vector(dx, dy, dz)

	def getr(self):
		if not self.r:
			self.r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
		return self.r
