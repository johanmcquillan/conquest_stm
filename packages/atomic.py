
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from packages.sph import sph
from packages.smartDict import SmartDict

class Radial(object):

	"""Stores the radial part of basis function and metadata,
	ie. quantum numbers (n and l) and zeta index.
	"""

	def __init__(self, zeta, n, l, radii, radialFuncValues, cutoff):
		"""Constructs radial part of a basis function.

		Args:
		    zeta (int): Indexes functions with the same n and l, but different cutoff
		    n (int): Principal quantum number
		    l (int): Orbital angular momentum quantum number
		    radii (float[]): List of radial distance values; Measured in Bohr radii
		    radialFuncValues (float[]): List of radial function values; Same length as r
		    cutoff (float): Range of radial function, beyond which value of R is 0.0
		"""
		self.zeta = zeta
		self.n = n
		self.l = l
		self.radii = radii
		self.radialFuncValues = radialFuncValues
		self.cutoff = cutoff

	def getValue(self, distance):
		"""Use linear interpolation to evaluate radial function at distance.
		
		Args:
		    distance (float): Distance from origin in Bohr radii
		
		Returns:
		    float: Value of radial part
		"""
		if distance > self.cutoff or distance > self.radii[-1]:
			value = 0.0
		else:
			# Find the first r value larger than distance
			i = 0
			while self.radii[i] < distance:
				i = i + 1

			i = min(range(0, len(self.radii)), key=lambda j: abs(self.radii[j] - distance))

			# Get nearest stored values
			x1 = self.radii[i-1]
			x2 = self.radii[i]
			y1 = self.radialFuncValues[i-1]
			y2 = self.radialFuncValues[i]

			# Calculate via interpolation
			value = y1 + (distance - x1) * (y2 - y1) / (x2 - x1)
		return value

class Ion(object):
	"""
	Stores data on an ion represented with a certain basis.

	Attributes:
		ionName (string): Name of ion (usually from name of .ion file)
		radials (SmartDict): Radial objects accessed by radials[zeta][n][l], where all indices are int

	"""
	def __init__(self, ionName, radialDict=SmartDict()):
		"""Constuct ion represented using given basis functions.
		
		Args:
		    ionName (string): Name of ion (usually from name of .ion file)
		    radialDict (SmartDict, optional): Radial objects, indexed by radialDict[zeta][n][l],
		    									where all indices are int. Default is empty, radials
		    									may be added after instantion
		"""
		self.ionName = ionName
		self.radials = radialDict # Radial objects; accessed by self.radials[zeta][n][l]

	def sortPAOs(self):
		"""Sort pseudo-atomic orbitals into order according to .dat files.
		
		Returns:
		    list: Ordered list of PAO data;
		    		Each element is a list containing [zeta, n, l, m] for the PAO
		"""
		sortedPAOs = []

		# Dict keys not necessarily in order
		# Order key list by lowest to highest for each index
		zetaList = sorted(self.radials.keys())
		for zeta in zetaList:
			nList = sorted(self.radials[zeta].keys())
			for n in nList:
				lList = sorted(self.radials[zeta][n])
				for l in lList:
					for m in range(-l, l+1):
						sortedPAOs.append([zeta, n, l, m])
		return sortedPAOs

	def hasRadial(self, zeta, n, l):
		"""Check if Ion has radial object for specified object without creating
		a dict.

		This is needed as SmartDict will automatically create an empty
		dict if given an invalid key.

		Args:
			zeta (int): Indexes functions with different cutoff but same n and l
		    n (int): Principal quantum number
		    l (int): Orbital angular momentum quantum number

		Returns:
			boolean: True if Radial is stored, false if not
		"""
		output = False
		if self.radials.has_key(zeta):
			if self.radials[zeta].has_key(n):
				if self.radials[zeta][n].has_key(l):
					return True
		return output

	def addRadial(self, radial):
		"""Adds radial to self.radials.

		Overwrites radial with same metadata (zeta, n, l).

		Args:
			radial (Radial): Radial object to add"""

		# Get metadata
		zeta = radial.zeta
		n = radial.n
		l = radial.l

		# Add Radial
		self.radials[zeta][n][l] = radial
		self.sortPAOs()

	def getRadial(self, zeta, n, l):
		"""Return Radial object for specified orbital.

		Encapsulation needed such that self.radials (SmartDict) does not create
		entry for invalid keys.

		Args:
			zeta (int): Indexes functions with different cutoff but same n and l
		    n (int): Principal quantum number
		    l (int): Orbital angular momentum quantum number

		Returns:
			Radial: Radial object for specified indices
		"""

		output = None
		if self.hasRadial(zeta, n, l):
			output = self.radials[zeta][n][l]
		return output

	def getRadialValue(self, zeta, n, l, r):
		"""Use linear interpolation to evaluate radial function at distance r.

		Encapsulation needed such that self.radials (SmartDict) does not create
		entry for invalid keys.

		Args:
			zeta (int): Indexes functions with different cutoff but same n and l
		    n (int): Principal quantum number
		    l (int): Orbital angular momentum quantum number
		    r (float): Radial distance from ion in Bohr radii

		Returns:
			float: Radial function evaluated at r
		"""

		output = None
		if self.hasRadial(zeta, n, l):
			output = self.radials[zeta][n][l].getValue(r)
		return output

	def getMaxCutoff(self):
		"""Return the maximum cutoff radius for all stored Radials as a float.

		Beyond the cutoff the radial part of the wavefunction is defined to be 0.		
		"""
		maxcut = 0.0
		for zeta in self.radials:
			for n in self.radials[zeta]:
				for l in self.radials[zeta][n]:
					if maxcut < self.radials[zeta][n][l].cutoff:
						maxcut = self.radials[zeta][n][l].cutoff
		return maxcut

class Atom(Ion):

	"""Stores information on an atom, extending Ion to include atomic position and basis coefficients

	Attributes:
		ionName (string): Name of ion (usually from name of .ion file)
		radials (SmartDict): Radial objects accessed by radials[zeta][n][l], where all indices are int
		x (float): Cartesian x-coordinate of atom
		y (float): Cartesian y-coordinate of atom
		z (float): Cartesian z-coordinate of atom
		bands (SmartDict): Nested dict of complex basis coefficients;
							Accessed by bands[bandEnergy][zeta][n][l][m]
	"""

	def __init__(self, ionName, x, y, z, radials=SmartDict()):
		"""Constructor for atom.

		Args:
			ionName (string): Name of ion (usually from name of .ion file)
		    radialDict (SmartDict, optional): Radial objects, indexed by radialDict[zeta][n][l],
		    									where all indices are int. Default is empty, radials
		    									may be added after instantiation
		    x (float): Atomic x coordinate in simulation cell in Bohr radii
		    y (float): Atomic y coordinate in simulation cell in Bohr radii
		    z (float): Atomic z coordinate in simulation cell in Bohr radii
		"""
		Ion.__init__(self, ionName, radials)
		self.x = x
		self.y = y
		self.z = z
		self.bands = SmartDict()

	def setIon(self, I):
		"""Copy all attributes from an Ion to this Atom.

		Args:
			I (Ion): Ion object from which to copy attributes
		"""
		self.ionName = I.ionName
		self.radials = I.radials
		self.sortPAOs()

	def hasCoeff(self, E, zeta, n, l, m):
		"""Check if atom stores coefficient for given orbital.

		Args:
			zeta (int): Indexes functions with different cutoff but same n and l
		    n (int): Principal quantum number
		    l (int): Orbital angular momentum quantum number
		    m (int): Azimuthal orbital angular momentum quantum number

		Returns:
			boolean: True if coefficient is stored, false if not
		"""

		output = False
		if self.bands.has_key(E):
			if self.bands[E].has_key(zeta):
				if self.bands[E][zeta].has_key(n):
					if self.bands[E][zeta][n].has_key(l):
						if self.bands[E][zeta][n][l].has_key(m):
							output = True
		return output

	def addCoeff(self, E, PAO, coeff):
		"""Add a complex coefficient to self.bands.

		Args:
			PAO (int): index of PAO as given in .dat file
			coeff (coomplex): Coefficient of PAO
		"""

		PAOdata = self.sortPAOs()[PAO - 1]
		zeta = PAOdata[0]
		n = PAOdata[1]
		l = PAOdata[2]
		m = PAOdata[3]

		self.bands[E][zeta][n][l][m] = coeff

	def getCoeff(self, E, zeta, n, l, m):
		"""Return complex coefficient for given orbital.

		Args:
			zeta (int): Indexes functions with different cutoff but same n and l
		    n (int): Principal quantum number
		    l (int): Orbital angular momentum quantum number
		    m (int): Azimuthal orbital angular momentum quantum number

		Returns:
			complex: Coefficient for given orbital
		"""
		output = None
		if self.hasCoeff(E, zeta, n, l, m):
			output = self.bands[E][zeta][n][l][m]
		return output

	def getPsi(self, E, x, y, z):
		"""Return complex wavefunction at (x,y,z) due to this atom only.

		Args:
			x (float): Cartesian x coordinate at which to evaluate wavefunction
			y (float): Cartesian y coordinate at which to evaluate wavefunction
			z (float): Cartesian z coordinate at which to evaluate wavefunction

		Returns:
			complex: Complex wavefunction evaluated at x, y, z
		"""

		# Get coords relative to atom
		relx = x - self.x
		rely = y - self.y
		relz = z - self.z
		# Get relative distance to atom
		r = np.sqrt(relx**2 + rely**2 + relz**2)

		psi = complex(0.0, 0.0)

		# If r is beyond atoms range, return 0.0+0.0j
		if r <= self.getMaxCutoff(): # Loop over all radial functions
			for zeta in self.radials:
				for n in self.radials[zeta]:
					for l in self.radials[zeta][n]:
						# Evaluate R(r) using linear interpolation
						R = self.getRadialValue(zeta, n, l, r)

						for m in range(-l, l+1):
							# Calculate spherical harmonic
							Y = sph(l, m, relx, rely, relz)

							# Get coefficient of basis functoin
							coeff = self.bands[E][zeta][n][l][m]
							# Calculate and add contribution of basis function
							psiReal = R*Y*coeff.real
							psiImag = R*Y*coeff.imag
							psi += complex(psiReal, psiImag)
		return psi
