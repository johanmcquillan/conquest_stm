import numpy as np
from scipy.interpolate import interp1d

from sph import sph
from smartDict import SmartDict


class Radial(object):
	"""Stores the radial part of basis function and metadata,
	ie. quantum numbers (n and l) and zeta index.

	All lengths measured in Bohr radii (a0).
	"""

	def __init__(self, n, l, zeta, radii, radialFuncValues, cutoff):
		"""Constructs radial part of a basis function.

		Args:
			zeta (int): Indexes functions with the same n and l, but different cutoff
			n (int): Principal quantum number
			l (int): Orbital angular momentum quantum number
			radii (float[]): List of radial distance values
			radialFuncValues (float[]): List of radial function values for each value of radii
			cutoff (float): Range of radial function, beyond which value of R is 0.0
		"""
		self.zeta = zeta
		self.n = n
		self.l = l
		self.radii = radii
		self.radialFuncValues = radialFuncValues
		self.radialInterpolator = interp1d(radii, radialFuncValues, kind='cubic')
		self.cutoff = cutoff

	def getValueLinear(self, distance):
		"""Use linear interpolation to evaluate radial function at distance.
		
		Args:
			distance (float): Distance from origin

		Returns:
			float: Value of radial part
		"""
		if distance > self.cutoff:
			value = 0.0
		else:
			# Find the first r value larger than distance
			i = 0
			while self.radii[i] < distance:
				i += 1

			# Get nearest stored values
			x1 = self.radii[i - 1]
			x2 = self.radii[i]
			y1 = self.radialFuncValues[i - 1]
			y2 = self.radialFuncValues[i]

			# Calculate via interpolation
			value = y1 + (distance - x1) * (y2 - y1) / (x2 - x1)
		return value

	def getValueCubic(self, distance):
		"""Use cubic interpolation to evaluate radial function at distance.
		
		Args:
			distance (float): Distance from origin
		
		Returns:
			float: Value of radial part
		"""

		if distance > self.cutoff:
			value = 0.0
		else:
			value = self.radialInterpolator(distance)
		return value


class Ion(object):
	"""
	Stores data on an ion represented with a certain basis.

	All lengths measured in Bohr radii (a0).

	Attributes:
		ionName (string): Name of ion (usually from name of .ion file)
		radials (SmartDict): Radial objects accessed by radials[zeta][n][l], where all indices are int

	"""

	def __init__(self, ionName, radialDict=SmartDict()):
		"""Construct ion represented using given basis functions.
		
		Args:
			ionName (string): Name of ion (usually from name of .ion file)
			radialDict (SmartDict, optional): Radial objects, indexed by radialDict[zeta][n][l],
												where all indices are int. Default is empty, radials
												may be added after instantiation
		"""
		self.ionName = ionName
		self.radials = radialDict  # Radial objects; accessed by self.radials[zeta][n][l]

	def sortPAOs(self):
		"""Sort pseudo-atomic orbitals into order according to .dat files.
		
		Returns:
			list: Ordered list of PAO data;
					Each element is a list containing [l, zeta, m] for the PAO
		"""
		sortedPAOs = []

		# Dict keys not necessarily in order
		# Order key list by lowest to highest for each index
		lList = sorted(self.radials.keys())
		for l in lList:
			zetaList = sorted(self.radials[l])
			for zeta in zetaList:
				for m in range(-l, l + 1):
					sortedPAOs.append([l, zeta, m])
		return sortedPAOs

	def hasRadial(self, l, zeta):
		"""Check if Ion has radial object for specified object without creating
		a dict.

		Encapsulation required due to autovivification of SmartDict.

		Args:
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoffs

		Returns:
			boolean: True if Radial is stored, false if not
		"""
		output = False
		if l in self.radials:
			if zeta in self.radials[l]:
				return True
		return output

	def addRadial(self, radial):
		"""Adds radial to self.radials.

		Overwrites radial with same metadata (l, zeta).

		Args:
			radial (Radial): Radial object to add
		"""

		# Get metadata
		zeta = radial.zeta
		l = radial.l

		# Add Radial
		self.radials[l][zeta] = radial
		self.sortPAOs()

	def getRadial(self, l, zeta):
		"""Return Radial object for specified orbital.

		Encapsulation needed such that self.radials (SmartDict) does not create
		entry for invalid keys.

		Args:
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoffs

		Returns:
			Radial: Radial object for specified indices
		"""

		output = None
		if self.hasRadial(l, zeta):
			output = self.radials[l][zeta]
		return output

	def getRadialValue(self, l, zeta, r, interpolation='cubic'):
		"""Use linear interpolation to evaluate radial function at distance r.

		Encapsulation required due to autovivification of SmartDict.

		Args:
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoffs
			r (float): Radial distance from ion
			interpolation (string, opt.): Method of interpolation; possible arguments are 'cubic' (default) and 'linear'

		Returns:
			float: Radial function evaluated at r
		"""

		output = None
		if self.hasRadial(l, zeta):
			if interpolation == 'cubic':
				output = self.radials[l][zeta].getValueCubic(r)
			elif interpolation == 'linear':
				output = self.radials[l][zeta].getValueLinear(r)
		return output

	def getMaxCutoff(self):
		"""Return the maximum cutoff radius for all stored Radials as a float.

		Beyond the cutoff the radial part of the wavefunction is defined to be 0.		
		"""
		maxCut = 0.0
		for l in self.radials:
			for zeta in self.radials[l]:
				if maxCut < self.radials[l][zeta].cutoff:
					maxCut = self.radials[l][zeta].cutoff
		return maxCut


class Atom(Ion):
	"""Stores information on an atom, extending Ion to include atomic position and basis coefficients

	All lengths measured in Bohr radii (a0).
	All energies measured in electron volts (eV).

	Attributes:
		ionName (string): Name of ion (usually from name of .ion file)
		radials (SmartDict): Radial objects accessed by radials[l][zeta], where all indices are int
		x (float): Cartesian x-coordinate of atom
		y (float): Cartesian y-coordinate of atom
		z (float): Cartesian z-coordinate of atom
		bands (SmartDict): Nested dict of complex basis coefficients;
							Accessed by bands[bandEnergy][l][zeta][m]
	"""

	def __init__(self, ionName, x, y, z, radials=SmartDict()):
		"""Constructor for atom.

		Args:
			ionName (string): Name of ion (usually from name of .ion file)
			radials (SmartDict, optional): Radial objects, indexed by radialDict[zeta][n][l],
											where all indices are int. Default is empty, radials
											may be added after instantiation
			x (float): Atomic x coordinate in simulation cell
			y (float): Atomic y coordinate in simulation cell
			z (float): Atomic z coordinate in simulation cell
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

	def hasCoefficient(self, Kx, Ky, Kz, E, l, zeta, m):
		"""Check if atom stores coefficient for given orbital.

		Encapsulation required due to autovivification of SmartDict.

		Args:
			Kx (float): K-point x coordinate
			Ky (float): K-point y coordinate
			Kz (float): K-point z coordinate
			E (float): Band energy
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoff but same n and l
			m (int): Azimuthal orbital angular momentum quantum number

		Returns:
			boolean: True if coefficient is stored, false if not
		"""

		output = False
		if Kx in self.bands:
			if Ky in self.bands[Kx]:
				if Kz in self.bands[Kx][Ky]:
					if E in self.bands[Kx][Ky][Kz]:
						if l in self.bands[Kx][Ky][Kz][E]:
							if zeta in self.bands[Kx][Ky][Kz][E][l]:
								if m in self.bands[Kx][Ky][Kz][E][l][zeta]:
									output = True
		return output

	def addCoefficient(self, Kx, Ky, Kz, E, PAO, coefficient):
		"""Add a complex coefficient to self.bands.

		Args:
			Kx (float): K-point x coordinate
			Ky (float): K-point y coordinate
			Kz (float): K-point z coordinate
			E (float): Band energy
			PAO (int): index of PAO as given in .dat file
			coefficient (complex): Coefficient of PAO
		"""

		PAOdata = self.sortPAOs()[PAO - 1]
		l = PAOdata[0]
		zeta = PAOdata[1]
		m = PAOdata[2]

		self.bands[Kx][Ky][Kz][E][l][zeta][m] = coefficient

	def getCoefficient(self, Kx, Ky, Kz, E, l, zeta, m):
		"""Return complex coefficient for given orbital.

		Args:
			Kx (float): K-point x coordinate
			Ky (float): K-point y coordinate
			Kz (float): K-point z coordinate
			E (float): Band energy
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoff but same n and l
			m (int): Azimuthal orbital angular momentum quantum number

		Returns:
			complex: Coefficient for given orbital
		"""
		output = 0.0
		if self.hasCoefficient(Kx, Ky, Kz, E, l, zeta, m):
			output = self.bands[Kx][Ky][Kz][E][l][zeta][m]
		return output

	def getTotalKPoints(self):
		"""Count total number of k-points"""
		totalKPoints = 0
		for Kx in self.bands:
			for Ky in self.bands[Kx]:
				for Kz in self.bands[Kx][Ky]:
					totalKPoints += 1
		return totalKPoints

	def getPsi(self, Kx, Ky, Kz, E, x, y, z, interpolation='cubic'):
		"""Evaluate wavefunction contribution from this atom.

		Args:
			Kx (float): K-point x coordinate
			Ky (float): K-point y coordinate
			Kz (float): K-point z coordinate
			E (float): Band energy
			x (float): Cartesian x coordinate
			y (float): Cartesian y coordinate
			z (float): Cartesian z coordinate
			interpolation (string, opt.): Method of interpolation; possible arguments are 'cubic' (default) and 'linear'

		Returns:
			complex: Wavefunction value
			"""
		psi = complex(0.0, 0.0)

		# Get relative displacement to atom
		relx = x - self.x
		rely = y - self.y
		relz = z - self.z
		# Get relative distance to atom
		r = np.sqrt(relx**2 + rely**2 + relz**2)
		# If r is beyond atoms range, return 0j
		if r <= self.getMaxCutoff():  # Loop over all radial functions
			for l in self.radials:
				for zeta in self.radials[l]:
					# Evaluate radial part
					R = self.getRadialValue(l, zeta, r, interpolation=interpolation)
					for m in range(-l, l + 1):
						# Evaluate spherical harmonic
						Y = sph(l, m, relx, rely, relz)
						# Get coefficient of basis function
						coefficient = self.getCoefficient(Kx, Ky, Kz, E, l, zeta, m)
						# Calculate and add contribution of basis function
						psi += R * Y * coefficient
		return psi

	def applyFactor(self, factor, Kx, Ky, Kz, E):
		"""Apply a normalisation factor to all coefficients for a given energy and k-point.

		Args:
			factor (float): Normalisation factor to be applied
			Kx (float): K-point x coordinate
			Ky (float): K-point y coordinate
			Kz (float): K-point z coordinate
			E (float): Energy of band to which to apply factor
		"""
		# Loop over all coefficients
		for l in self.bands[Kx][Ky][Kz][E]:
			for zeta in self.bands[Kx][Ky][Kz][E][l]:
				for m in self.bands[Kx][Ky][Kz][E][l][zeta]:
					# Apply factor
					self.bands[Kx][Ky][Kz][E][l][zeta][m] *= factor
