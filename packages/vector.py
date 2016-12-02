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
		self.magnitude = None

	def __str__(self):
		return '('+str(self.x)+', '+str(self.y)+', '+str(self.z)+')'

	def __hash__(self):
		return hash((self.x, self.y, self.z))

	def __eq__(self, other):
		return (self.x, self.y, self.z) == (other.x, other.y, other.z)

	def __ne__(self, other):
		return not self == other

	def subtract(self, V):
		"""Subtract vector from this vector"""
		dx = self.x - V.x
		dy = self.y - V.y
		dz = self.z - V.z
		return Vector(dx, dy, dz)

	def getMagnitude(self):
		"""Get vector magnitude"""
		if self == self.zero():
			self.magnitude = 0
		elif not self.magnitude:
			self.magnitude = np.sqrt(self.x**2 + self.y**2 + self.z**2)
		return self.magnitude

	@staticmethod
	def zero():
		"""Zero vector"""
		return Vector(0.0, 0.0, 0.0)
