
import packages.atomic

import numpy as np

class Cell(object):

	"""Simulation cell which holds Atom objects in a 3D mesh."""

	def __init__(self, name, xLength, yLength, zLength, gridSpacing=0.1):
		"""Constructs 3D cell with given dimensional.

		Args:
			xLength (float): Length of cell along x in Bohr radii
			yLength (float): Length of cell along y in Bohr radii
			zLength (float): Length of cell along z in Bohr radii
			gridSpacing (float optional): Resolution of mesh points
		"""
		self.name = name
		self.xLength = xLength
		self.yLength = yLength
		self.zLength = zLength
		self.gridSpacing = gridSpacing

		self.xPoints = int(xLength / gridSpacing)
		self.yPoints = int(yLength / gridSpacing)
		self.zPoints = int(zLength / gridSpacing)

		self.xMesh, self.yMesh, self.zMesh = np.mgrid[0:xLength:gridSpacing,
		                                              0:yLength:gridSpacing,
		                                              0:zLength:gridSpacing]
		self.psi = np.empty_like(self.xMesh, dtype=complex)
		self.atoms = {}
		self.bands = []

	def addAtom(self, atom, atomKey):
		"""Add atom to self.atoms, indexed by atomKey"""

		# Reassign atom coordinates to nearest mesh points
		atom.x = self.gridSpacing*min(range(0, self.xPoints),
			                             key=lambda i: abs(self.xMesh[i, 0, 0]-atom.x))
		atom.y = self.gridSpacing*min(range(0, self.yPoints),
			                             key=lambda i: abs(self.yMesh[0, i, 0]-atom.y))
		atom.z = self.gridSpacing*min(range(0, self.zPoints),
			                             key=lambda i: abs(self.zMesh[0, 0, i]-atom.z))

		# Add to dict
		self.atoms[atomKey] = atom
		for band in atom.bands:
			if band not in self.bands:
				self.bands.append(band)
		self.bands = sorted(self.bands)

	def setPsi(self, band=0, debug=False):
		"""Calculate complex wavefunction at all points in 3D mesh and
		assign it to self.psi"""

		E=self.bands[band]

		# Get 3D mesh with (0+0j) at all points
		wavefunc = np.empty_like(self.xMesh, dtype=complex)
		# Iterate over all atoms stored in this cell
		for atomKey in self.atoms:
			atom = self.atoms[atomKey]
			# Iterate over all mesh points
			for i in range(0, self.xPoints):
				for j in range(0, self.yPoints):
					for k in range(0, self.zPoints):
						# Get mesh coordinates
						x = self.xMesh[i, j, k]
						y = self.yMesh[i, j, k]
						z = self.zMesh[i, j, k]

						# # Add contribution from this atom to this mesh point
						wavefunc[i, j, k] = wavefunc[i, j, k] + atom.getPsi(E, x, y, z)
			if debug:
				print 'Band Energy = '+str(E)+'; Calculated Psi for Atom '+str(atomKey)
		self.psi = wavefunc

	def givePsi(self, x, y, z, band=0, debug=False):

		psi = complex(0.0, 0.0)
		E = self.bands[band]
		for atomKey in self.atoms:
			atom = self.atoms[atomKey]
			if atom.bands.has_key(E):
				psi = psi + atom.getPsi(E, x, y, z)
			else:
				print "Atom "+str(a)
		return psi



