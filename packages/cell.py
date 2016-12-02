import numpy as np
from smartDict import SmartDict
from vector import Vector

BOLTZMANN = 8.6173303E-5  # Boltzmann's Constant in eV/K


class Cell(object):
	"""Simulation cell which holds Atom objects in a 3D mesh.

	All lengths measured in Bohr radii (a0).
	All energies measured in Hartrees (Ha).
	
	Attributes:
		name (string): Name of simulation; used for plot titles
			fermiLevel (float): Fermi Level of simulation
			xLength (float): Length of cell along x
			yLength (float): Length of cell along y
			zLength (float): Length of cell along z
			gridSpacing (float): Resolution of mesh points
			xPoints (float): Number of points on x mesh
			yPoints (float): Number of points on y mesh
			zPoints (float): Number of points on z mesh
			xMesh (3D mgrid): Mesh of x values
			yMesh (3D mgrid): Mesh of y values
			zMesh (3D mgrid): Mesh of z values
			atoms (int : Atom): Atom objects of simulation indexed by atom number
			basisPoints (SmartDict): Basis function values, indexed by [atomKey][x][y][z][l][zeta][m]
			bands (SmartDict): Energies of bands indexed by k-point (eg. [Kx][Ky][Kz])
	"""

	def __init__(self, name, fermiLevel, electrons, xLength, yLength, zLength, gridSpacing=0.1):
		"""Constructs 3D cell with given dimensional.

		All lengths measured in Bohr radii (a0);
		All energies measured in Hartrees (Ha).

		Args:
			name (string): Name of simulation; used for plot titles
			fermiLevel (float): Fermi Level of simulation
			xLength (float): Length of cell along x
			yLength (float): Length of cell along y
			zLength (float): Length of cell along z
			gridSpacing (float, opt.): Resolution of mesh points
		"""

		self.name = name
		self.fermiLevel = fermiLevel
		self.electrons = electrons
		self.xLength = xLength
		self.yLength = yLength
		self.zLength = zLength
		self.gridSpacing = gridSpacing

		# Calculator number of points
		self.xPoints = int(xLength / gridSpacing)
		self.yPoints = int(yLength / gridSpacing)
		self.zPoints = int(zLength / gridSpacing)

		# Form Cartesian meshes
		self.xMesh, self.yMesh, self.zMesh = np.mgrid[0: xLength: gridSpacing, 0: yLength: gridSpacing, 0: zLength: gridSpacing]

		# Initialise atoms and bands
		self.atoms = {}
		self.basisPoints = SmartDict()
		self.bands = SmartDict()

	def hasKPoint(self, Kx, Ky, Kz):
		"""Check if cell stores specified k-point.

		Encapsulation required due to autovivification of SmartDict.

		Args:
			Kx (float): K-point x coordinate
			Ky (float): K-point y coordinate
			Kz (float): K-point z coordinate
		"""
		output = False
		if Kx in self.bands:
			if Ky in self.bands[Kx]:
				if Kz in self.bands[Kx][Ky]:
					output = True
		return output

	def hasBand(self, Kx, Ky, Kz, E):
		"""Check if cell stores specified band.

		Encapsulation required due to autovivification of SmartDict.

		Args:
			Kx (float): K-point x coordinate
			Ky (float): K-point y coordinate
			Kz (float): K-point z coordinate
			E (float): Band energy
		"""
		output = False
		if self.hasKPoint(Kx, Ky, Kz):
			if E in self.bands[Kx][Ky][Kz]:
				output = True
		return output

	def setBasisPoint(self, atomKey, position, interpolation='cubic'):
		"""Calculate and save basis function values for all points within cutoff region.

		Args:
			atomKey (int): Atom number, as given in Conquest_out
			position (Vector): 3D Cartesian real space vector
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
		"""
		atom = self.atoms[atomKey]
		Y = 0.0
		# Check if atom is in range
		if self.atoms[atomKey].withinCutoff(position):
			for l in atom.radials:
				for zeta in atom.radials[l]:
					R = atom.getRadialValueRelative(l, zeta, position, interpolation=interpolation)
					for m in range(-l, l+1):
						if R != 0:
							Y = atom.getSPHrelative(l, m, position)
						self.basisPoints[atomKey][position.x][position.y][position.z][l][zeta][m] = R*Y
		else:
			self.basisPoints[atomKey][position.x][position.y][position.z] = SmartDict()

	def hasBasisPoint(self, atomKey, position):
		"""Check if basis point has already been calculated

		Args:
			atomKey (int): Atom number, as given in Conquest_out
			position (Vector): 3D Cartesian real space vector
		"""
		output = False
		if atomKey in self.basisPoints:
			if position.x in self.basisPoints[atomKey]:
				if position.y in self.basisPoints[atomKey][position.x]:
					if position.z in self.basisPoints[atomKey][position.x][position.y]:
						output = True
		return output

	def addAtom(self, atom, atomKey):
		"""Add atom to self.atoms, indexed by atomKey

		Args:
			atom (Atom): Atom object
			atomKey (int): Atom number, as given in Conquest_out
		"""
		# Add to dict
		self.atoms[atomKey] = atom

		# Loop over atoms k-points
		for Kx in atom.bands:
			for Ky in atom.bands[Kx]:
				for Kz in atom.bands[Kx][Ky]:
					# If cell does not have k-point, create empty band energy list
					if not self.hasKPoint(Kx, Ky, Kz):
						self.bands[Kx][Ky][Kz] = []
					# Add band energies to k-point
					for E in atom.bands[Kx][Ky][Kz]:
						self.bands[Kx][Ky][Kz].append(E)
					# Sort energy list
					self.bands[Kx][Ky][Kz] = sorted(self.bands[Kx][Ky][Kz])

	def getGammaEnergies(self):
		"""Return list of energies at gamma-point"""
		return sorted(self.bands[0.0][0.0][0.0])

	def getTotalKPoints(self):
		"""Count total number of k-points"""
		totalKPoints = 0
		for Kx in self.bands:
			for Ky in self.bands[Kx]:
				for Kz in self.bands[Kx][Ky]:
					totalKPoints += 1
		return totalKPoints

	def fermiDirac(self, energy, temperature):
		"""Calculate Fermi-Dirac distribution value.

		Args:
			energy (float): Energy in eV
			temperature (float): Absolute temperature in K

		Returns:
			float: Occupation value
		"""
		f = 1.0 / (np.exp((energy - self.fermiLevel) / (BOLTZMANN * temperature)) + 1)
		return f

	def getPsi(self, Kx, Ky, Kz, E, position, interpolation='cubic'):
		"""Evaluate wavefunction at specific position, k-point, and energy.

		Args:
			Kx (float): K-point x coordinate
			Ky (float): K-point y coordinate
			Kz (float): K-point z coordinate
			E (float): Band energy
			position (Vector): 3D Cartesian real space vector
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'

		Returns:
			complex: Wavefunction value
		"""
		psi = complex(0.0, 0.0)
		# Loop over atoms
		for atomKey in self.atoms:
			# Check is basis function values have been calculated for this atom and point
			if not self.hasBasisPoint(atomKey, position):
				# Get and store basis function values
				# Store values so they only need to be calculated once
				self.setBasisPoint(atomKey, position, interpolation=interpolation)
			# Use basis function value to calculate psi
			BP = self.basisPoints[atomKey][position.x][position.y][position.z]
			psi += self.atoms[atomKey].getPsi(Kx, Ky, Kz, E, position, basisPoint=BP, local=True)
		return psi

	def getPsiGamma(self, E, position):
		"""Evaluate wavefunction at specific position and energy using only gamma-point.

		Args:
			E (float): Band energy
			position (Vector): 3D Cartesian real space vector

		Returns:
			complex: Wavefunction value
			"""
		return self.getPsi(0.0, 0.0, 0.0, E, position)

	def getLDoS(self, Emin, Emax, T, position, interpolation='cubic', debug=False):
		"""Evaluate local density of states (LDoS) within energy range at specific point.

		Args:
			Emin (float): Minimum energy
			Emax (float): Maximum energy
			T (float): Absolute temperature in K
			position (Vector): 3D Cartesian real space vector
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
			debug (bool, opt.): If true, print extra information during runtime

		Returns:
			float: LDoS value
		"""
		I = 0.0
		totalK = self.getTotalKPoints()
		w = 1.0/totalK  # K-point weighting

		# Loop over all k-points
		for Kx in self.bands:
			for Ky in self.bands[Kx]:
				for Kz in self.bands[Kx][Ky]:
					# Loop over bands in energy range
					for E in self.bands[Kx][Ky][Kz]:
						if Emin < E < Emax:
							# Calculate LDoS
							psi = self.getPsi(Kx, Ky, Kz, E, position, interpolation=interpolation)
							I += w * self.fermiDirac(E, T) * (abs(psi))**2
		if debug:
			print "Finished LDoS = "+str(I)+", at r = "+str(position)
		return I
