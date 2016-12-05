import numpy as np


class Vector(object):
	"""3D cartesian vector class.

	Attributes:
		x (float): Component value in x direction
		y (float): Component value in y direction
		z (float): Component value in z direction
		magnitude (float): Vector magnitude. Initially null - only calculates if needed, saves computational cost
	"""

	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
		self.magnitude = None

	def __str__(self):
		"""Return string of components"""
		return '('+str(self.x)+', '+str(self.y)+', '+str(self.z)+')'

	def __hash__(self):
		"""Return hash of tuple of components"""
		return hash((self.x, self.y, self.z))

	def __eq__(self, other):
		"""Check if each component of vectors are equal"""
		return (self.x, self.y, self.z) == (other.x, other.y, other.z)

	def __ne__(self, other):
		"""Inverse of self -- other"""
		return not self == other

	def __abs__(self):
		"""Return magnitude of vector.

		self.magnitude is not initially calculated to save time.
		After initial calculation, it is saved for future reference.
		"""
		if self == self.zero():
			self.magnitude = 0
		elif not self.magnitude:
			self.magnitude = np.sqrt(self.x**2 + self.y**2 + self.z**2)
		return self.magnitude

	def __neg__(self):
		"""Return negative of vector"""
		return Vector(-self.x, -self.y, -self.z)

	def __mul__(self, other):
		"""Multiply vector by number"""
		return Vector(self.x * other, self.y * other, self.z * other)

	def __rmul__(self, other):
		"""Multiply vector by number"""
		return self * other

	def __div__(self, other):
		"""Divide vector by number"""
		return self * (1.0/other)

	def __add__(self, other):
		"""Add two vectors"""
		return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

	def __sub__(self, other):
		"""Subtract two vectors"""
		return self + (-other)

	def project_x(self):
		"""Project vector onto x-axis"""
		return Vector(self.x, 0, 0)

	def project_y(self):
		"""Project vector onto y-axis"""
		return Vector(0, self.y, 0)

	def project_z(self):
		"""Project vector onto z-axis"""
		return Vector(0, 0, self.z)

	def dot(self, other):
		"""Return dot product of this vector and an other vector"""
		return self.x*other.x + self.y*other.y + self.z*other.z

	def cross(self, other):
		"""Return cross product of this vector and an other vector"""
		x = self.y*other.z - self.z*other.y
		y = self.z*other.x - self.x*other.z
		z = self.x*other.y - self.y*other.x
		return Vector(x, y, z)

	def normalise(self):
		"""Return vector normalised to unit magnitude"""
		magnitude = abs(self)
		if magnitude != 0:
			return self / magnitude
		else:
			return self

	def is_positive(self):
		"""Return true if all components are greater or equal to 0"""
		if self.x >= 0 and self.y >= 0 and self.z >= 0:
			return True
		else:
			return False

	@staticmethod
	def zero():
		"""Zero vector"""
		return Vector(0.0, 0.0, 0.0)
