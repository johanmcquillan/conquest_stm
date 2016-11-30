import numpy as np
from smartDict import SmartDict

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
			atoms (int : Atom): Atom objects of simulation
								indexed by atom number
			bands ([float]): Energies of bands in ascending order
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
		self.bands = SmartDict()

	def hasKPoint(self, Kx, Ky, Kz):
		"""Check if cell stores specified k-point.

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

	def addAtom(self, atom, atomKey):
		"""Add atom to self.atoms, indexed by atomKey

		Args:
			atom (Atom): Atom object to add to simulation
			atomKey (int): Atom number, as given in Conquest_out
		"""

		# Reassign atom coordinates to nearest mesh points
		atom.x = self.gridSpacing * min(range(0, self.xPoints), key=lambda i: abs(self.xMesh[i, 0, 0] - atom.x))
		atom.y = self.gridSpacing * min(range(0, self.yPoints), key=lambda i: abs(self.yMesh[0, i, 0] - atom.y))
		atom.z = self.gridSpacing * min(range(0, self.zPoints), key=lambda i: abs(self.zMesh[0, 0, i] - atom.z))

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

	def getWavefunction(self, Emin, Emax, x, y, z):
		"""Evaluate complex wavefunction at specific point over specified energy range.

		Args:
			Emin (float): Minimum of energy range
			Emax (float): Maximum of energy range
			x (float): Cartesian x coordinate
			y (float): Cartesian y coordinate
			z (float): Cartesian z coordinate

		Returns:
			complex: Wavefunction value
		"""

		# Initialise psi and get band energy
		psi = complex(0.0, 0.0)

		# Get wavefunction contribution from each atom and add to psi
		for atomKey in self.atoms:
			atom = self.atoms[atomKey]
			psi = psi + atom.getPsi(Emin, Emax, x, y, z)
		return psi

	def getTotalCharge(self, Emin, Emax):
		"""Perform volume integration charge density of band over cell.

		Args:
			Emin (float): Minimum of energy range
			Emax (float): Maximum of energy range

		Returns:
			float: Total charge
		"""
		totalCharge = 0.0

		# Iterate over all mesh points
		for i in range(0, self.xPoints):
			for j in range(0, self.yPoints):
				for k in range(0, self.zPoints):
					# Get mesh coordinates
					x = self.xMesh[i, j, k]
					y = self.yMesh[i, j, k]
					z = self.zMesh[i, j, k]
					psi = 0.0
					# Sum contributions to psi from all atoms
					psi += self.getWavefunction(Emin, Emax, x, y, z)
					# Add to charge
					totalCharge += abs(psi)**2 * self.gridSpacing**3
		return totalCharge

	def normaliseBand(self, Emin, Emax, debug=False):
		"""Normalise coefficients of a band to give total charge of unity.

		Args:
			Emin (float): Minimum of energy range
			Emax (float): Maximum of energy range
			debug (bool, opt.): If true, print extra information during runtime
		"""
		totalCharge = self.getTotalCharge(Emin, Emax)

		# Check if difference between actual and calculated charge is large than tolerance
		if abs(1.0 - totalCharge) > 0.001:
			if debug:
				print "Total Electron Charge Unnormalised = " + str(totalCharge)

			# Apply normalisation factor to basis coefficients
			factor = np.sqrt(1.0 / totalCharge)
			for atomKey in self.atoms:
				self.atoms[atomKey].applyFactorForRange(factor, Emin, Emax)
		elif debug:
			print "Total Electron Charge Already Normalised"

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
			fermiLevel (float): Fermi Level in eV
			temperature (float): Absolute temperature in K

		returns:
			float: Occupation value
		"""
		f = 1.0 / (np.exp((energy - self.fermiLevel) / (BOLTZMANN * temperature)) + 1)
		return f

	def getCurrent(self, Emin, Emax, T, x, y, z, debug=False):

		I = 0.0
		w = 1.0/self.getTotalKPoints()

		for Kx in self.bands:
			for Ky in self.bands[Kx]:
				for Kz in self.bands[Kx][Ky]:
					for E in self.bands[Kx][Ky][Kz]:
						if Emin < E < Emax:
							psi = complex(0.0, 0.0)
							for atomKey in self.atoms:
								atom = self.atoms[atomKey]
								psi += atom.getPsi1(Kx, Ky, Kz, E, x, y, z)
							I += w * self.fermiDirac(E, T) * (abs(psi))**2
					if debug:
						print "Finished K-Point = "+str(Kx)+", "+str(Ky)+", "+str(Kz)
		print "Finished current I = "+str(I)+", at r = ("+str(x)+", "+str(y)+", "+str(z)+")"
		return I
