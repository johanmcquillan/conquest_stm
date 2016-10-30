
import numpy as np

class Vector:


	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def __str__(self):
		return '('+str(self.x)+', '+str(self.y)+', '+str(self.z)+')'

	def getComponants(self):
		return [self.x, self.y, self.z]


	def mag(self):
		r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
		return r

	def getUnitV(self):
		m = self.mag()
		x = self.x / m
		y = self.y / m
		z = self.z / m
		return Vector(x, y, z)

	def normalise(self):
		m = self.mag()
		self.x = self.x / m
		self.y = self.y / m
		self.z = self.z / m

