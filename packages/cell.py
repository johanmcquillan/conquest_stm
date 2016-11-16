
import packages.atomic

import numpy as np

class Cell(object):

	"""Simulation cell which holds Atom objects in a 3D mesh.
	
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

		All lengths measured in bohr radii (a0);
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

		# Calculater number of points
		self.xPoints = int(xLength / gridSpacing)
		self.yPoints = int(yLength / gridSpacing)
		self.zPoints = int(zLength / gridSpacing)

		# Form Cartesian meshes
		self.xMesh, self.yMesh, self.zMesh = np.mgrid[0:xLength:gridSpacing,
		                                              0:yLength:gridSpacing,
		                                              0:zLength:gridSpacing]

		# Initialise atoms and bands
		self.atoms = {}
		self.bands = []

	def addAtom(self, atom, atomKey):
		"""Add atom to self.atoms, indexed by atomKey

		Args:
			atom (Atom): Atom object to add to simulation
			atomKey (int): Atom number, as given in Conquest_out
		"""

		# Reassign atom coordinates to nearest mesh points
		atom.x = self.gridSpacing*min(range(0, self.xPoints),
			                             key=lambda i: abs(self.xMesh[i, 0, 0]-atom.x))
		atom.y = self.gridSpacing*min(range(0, self.yPoints),
			                             key=lambda i: abs(self.yMesh[0, i, 0]-atom.y))
		atom.z = self.gridSpacing*min(range(0, self.zPoints),
			                             key=lambda i: abs(self.zMesh[0, 0, i]-atom.z))

		# Add to dict
		self.atoms[atomKey] = atom

		# Add band energy if not already stored
		for band in atom.bands:
			if band not in self.bands:
				self.bands.append(band)
		self.bands = sorted(self.bands)

	def getPsiMesh(self, bandNumber=0, debug=False):
		"""Calculate complex wavefunction at all points in 3D mesh and
		assign it to self.psi

		band (int, opt.): Number of band at which to evaluate psi
		debug (bool., opt.): If true, print extra information during runtime
		
		"""

		bandEnergy=self.bands[bandNumber]

		# Get 3D mesh with (0+0j) at all points
		wavefunc = np.empty_like(self.xMesh, dtype=complex)
		# Iterate over all atoms stored in this cell
		for atomKey in self.atoms:
			atom = self.atoms[atomKey]
			if atom.bands.has_key(bandEnergy):
				# Iterate over all mesh points
				for i in range(0, self.xPoints):
					for j in range(0, self.yPoints):
						for k in range(0, self.zPoints):
							# Get mesh coordinates
							x = self.xMesh[i, j, k]
							y = self.yMesh[i, j, k]
							z = self.zMesh[i, j, k]

							# # Add contribution from this atom to this mesh point
							wavefunc[i, j, k] = wavefunc[i, j, k] + atom.getPsi(bandEnergy, x, y, z)
				if debug:
					print 'Band Energy = '+str(bandEnergy)+'; Calculated Psi for Atom '+str(atomKey)
			elif debug:
				print "Atom "+str(a)+" has no band "+str(bandNumber)+": "+str(bandEnergy)
		return wavefunc

	def givePsi(self, x, y, z, bandNumber=0, debug=False):
		"""Evaluate complex wavefunction at specific point.

		Args:
			x (float): Cartesian x coordinate
			y (float): Cartesian y coordinate
			z (float): Cartesian z coordinate
			bandNumber (int, opt.): Band number at which to evaluate
			debug (bool., opt.): If true, print extra information during runtime

		Returns:
			complex: Wavefunction value
		"""

		# Initialise psi and get band energy
		psi = complex(0.0, 0.0)
		bandEnergy = self.bands[bandNumber]

		# Get wavefunction contribution from each atom and add to psi
		for atomKey in self.atoms:
			atom = self.atoms[atomKey]
			if atom.bands.has_key(bandEnergy):
				psi = psi + atom.getPsi(bandEnergy, x, y, z)
			elif debug:
				print "Atom "+str(a)+" has no band "+str(bandNumber)+": "+str(bandEnergy)
		
		return psi

	def normaliseBand(self, bandNumber, debug=False):
		"""Normalise coefficients of a band to give a total electron charge of unity.

		Args:
			bandNumber (int): Number of band to normalise
			debug (bool, opt.): If true, print extra information during runtime
		"""

		totalCharge = 0.0
		bandEnergy=self.bands[bandNumber]

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
					for atomKey in self.atoms:
						atom = self.atoms[atomKey]
						if atom.bands.has_key(bandEnergy):
							psi += atom.getPsi(bandEnergy, x, y, z)
					# Add to charge
					totalCharge += abs(psi)**2 * self.gridSpacing**3
		
		# Check if difference between actual and calculated charge is large than tolerance
		if abs(self.electrons - totalCharge) > 0.001:
			if debug:
				print "Total Electron Charge Unnormalised = "+str(totalCharge)

			# Apply normalisation factor to basis coefficients
			factor =  np.sqrt(float(self.electrons) / totalCharge)
			for atomKey in self.atoms:
				self.atoms[atomKey].applyFactor(factor, bandEnergy)

		elif debug:
			print "Total Electron Charge Already Normalised"

	def combineBands(self, Erange, normalise=True, outputBandData=False, debug=False):
		"""Combine bands with energies within Erange.

		"""

		Emin = Erange[0]
		Emax = Erange[1]
		Eavg = (Emax - Emin) / 2.0

		bandsToCombine = []

		for i in range(len(self.bands)):
			if Emin < self.bands[i] and Emax > self.bands[i]:
				bandsToCombine.append(i)

		if len(bandsToCombine) > 2:
			for band in bandsToCombine:
				bandEnergy = self.bands[band]

				if normalise:
					self.normaliseBand(band, debug=debug)
					if debug:
						print 'Normalised band '+str(bandEnergy)+' Ha'

				for atomKey in self.atoms:
					self.atoms[atomKey].combineCoeffs(bandEnergy, Eavg)
				self.bands.remove(bandEnergy)

			self.bands.append(Eavg)
			self.bands = sorted(self.bands)
			i = 0
			newBandNumber = None
			bandFound = False
			print Eavg
			while not bandFound and i < len(self.bands):
				print self.bands[i]
				if self.bands[i] == Eavg:
					newBandNumber = i
					bandFound = True
				else:
					i += 1
			self.normaliseBand(newBandNumber)

		elif debug:
			if len(bandsToCombine) == 1:
				debugString = 'Only single band'
			else:
				debugString = 'No bands'
			print debugString+' within range '+str(Emin)+' Ha and '+str(Emax)+' Ha'

		if outputBandData:
			return newBandNumber, Eavg
