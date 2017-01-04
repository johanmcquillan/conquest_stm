import numpy as np


class Vector(object):
	"""3D real space Cartesian vector.

	Attributes:
		x (float): Component value in x direction
		y (float): Component value in y direction
		z (float): Component value in z direction
		magnitude (float): Vector magnitude; Initially null - only calculates if needed, saves computational cost
		hash (int): Hash value; Initially null - only calculates if needed, saves computational cost
	"""

	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
		self.magnitude = None
		self.hash = None

	def __str__(self):
		"""Return string of components"""
		return '('+str(self.x)+', '+str(self.y)+', '+str(self.z)+')'

	def __hash__(self):
		"""Return hash of tuple of components"""
		if self.hash is None:
			self.hash = hash((self.x, self.y, self.z))
		return self.hash

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
		elif self.magnitude is None:
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

	def reset_hash_mag(self):
		self.hash = None
		self.magnitude = None

	def constrain_vector_relative(self, cell_vector):
		"""Constrain vector to lie within simulation cell"""

		# Check if vector components are greater than half of cell sides
		# If greater, add or subtract cell length
		while self.x > cell_vector.x/2:
			self.x -= cell_vector.x
		while self.x <= -cell_vector.x/2:
			self.x += cell_vector.x

		while self.y > cell_vector.y/2:
			self.y -= cell_vector.y
		while self.y <= -cell_vector.y/2:
			self.y += cell_vector.y

		while self.z > cell_vector.z/2:
			self.z -= cell_vector.z
		while self.z <= -cell_vector.z/2:
			self.z += cell_vector.z

		self.reset_hash_mag()

	def constrain_vector_to_cell(self, cell_vector):
		"""Constrain vector to lie within simulation cell"""

		# Check if vector components are greater than half of cell sides
		# If greater, add or subtract cell length
		while self.x >= cell_vector.x:
			self.x -= cell_vector.x
		while self.x < 0:
			self.x += cell_vector.x

		while self.y >= cell_vector.y:
			self.y -= cell_vector.y
		while self.y < 0:
			self.y += cell_vector.y

		while self.z >= cell_vector.z:
			self.z -= cell_vector.z
		while self.z < 0:
			self.z += cell_vector.z

		self.reset_hash_mag()

	@staticmethod
	def zero():
		"""Zero vector"""
		return Vector(0.0, 0.0, 0.0)


class KVector(Vector):
	"""3D k-sapce Cartesian vector.

	Attributes:
		x (float): Component value in x direction
		y (float): Component value in y direction
		z (float): Component value in z direction
		magnitude (float): Vector magnitude. Initially null - only calculates if needed, saves computational cost
		hash (int): Hash value; Initially null - only calculates if needed, saves computational cost
		weight (float): K-point weighting
	"""

	def __init__(self, x, y, z, weight):
		Vector.__init__(self, x, y, z)
		self.weight = weight

	@staticmethod
	def gamma():
		"""Gamma point vector"""
		return KVector(0.0, 0.0, 0.0, 1.0)

