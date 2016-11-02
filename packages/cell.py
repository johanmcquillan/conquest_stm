
import atomic as atm

import numpy as np

class Cell(object):

	'''Simulation cell which holds Atom objects in a 3D mesh.'''

	def __init__(self, xLength, yLength, zLength, gridSpacing=0.1):
		self.xLength = xLength
		self.yLength = yLength
		self.zLength = zLength
		self.gridSpacing = gridSpacing

		self.xPoints = int(xLength / gridSpacing)
		self.yPoints = int(yLength / gridSpacing)
		self.zPoints = int(zLength / gridSpacing)

		self.xMesh, self.yMesh, self.zMesh = np.mgrid[0.0:xLength:gridSpacing, 
													  0.0:yLength:gridSpacing,
													  0.0:zLength:gridSpacing]
		self.psi = np.empty_like(self.xMesh)
		self.atoms = {}

	def addAtom(self, A, index):

		A.x = self.gridSpacing*min(range(0, self.xPoints), key=lambda i: abs(self.xMesh[i,0,0]-A.x))
		A.y = self.gridSpacing*min(range(0, self.yPoints), key=lambda i: abs(self.yMesh[0,i,0]-A.y))
		A.z = self.gridSpacing*min(range(0, self.zPoints), key=lambda i: abs(self.zMesh[0,0,i]-A.z))

		self.atoms[index] = A

	def setPsi(self):

		psi = np.empty_like(self.xMesh, dtype=np.complex64)

		for a in self.atoms.keys():
			for i in range(0, self.xPoints):
				for j in range(0, self.yPoints):
					for k in range(0, self.zPoints):
						x = self.xMesh[i,j,k]
						y = self.yMesh[i,j,k]
						z = self.xMesh[i,j,k]

						psi[i,j,k] += self.atoms[a].getPsi(x, y, z)
			print 'Calculated atom '+str(a)
		self.psi = psi


	def plot(self):
		return None





