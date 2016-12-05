import numpy as np
from scipy.interpolate import interp1d

from sph import sph
from smartDict import SmartDict


class Radial(object):
	"""Stores the radial part of basis function and metadata,
	ie. quantum numbers (n and l) and zeta index.

	All lengths measured in Bohr radii (a0).
	"""

	interpolators = ['linear', 'quadratic', 'cubic']

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
		self.interpolator_quadratic = interp1d(radii, radialFuncValues, kind='quadratic')
		self.interpolator_cubic = interp1d(radii, radialFuncValues, kind='cubic')
		self.cutoff = cutoff

	def interpolator_linear(self, distance):
		"""Use linear interpolation to calculate radial function value at distance"""
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
		return y1 + (distance - x1) * (y2 - y1) / (x2 - x1)

	def get_value(self, distance, interpolation='cubic'):
		"""Evaluate radial function value at distance using interpolation method.

		Args:
			distance (float): Distance to evaluate radial function
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'

		Returns:
			float: Radial function value
		"""
		if distance > self.cutoff:
			value = 0.0
		else:
			if interpolation == 'linear':
				value = self.interpolator_linear(distance)
			elif interpolation == 'quadratic':
				value = self.interpolator_quadratic(distance)
			elif interpolation == 'cubic':
				value = self.interpolator_cubic(distance)
			else:
				errorString = "Interpolation method must be either "
				for i in range(len(self.interpolators)):
					errorString += self.interpolators[i]
					if i < len(self.interpolators) - 1:
						errorString += ', '
				raise ValueError(errorString)
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

	def has_radial(self, l, zeta):
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

	def add_radial(self, radial):
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

	def get_radial(self, l, zeta):
		"""Return Radial object for specified orbital.

		Encapsulation needed such that self.radials (SmartDict) does not create
		entry for invalid keys.

		Args:
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoffs

		Returns:
			Radial: Radial object for specified indices
		"""

		if self.has_radial(l, zeta):
			return self.radials[l][zeta]
		else:
			raise ValueError("No radial data for "+self.ionName+", l="+str(l)+", zeta="+str(zeta))

	def get_radial_value(self, l, zeta, distance, interpolation='cubic'):
		"""Use linear interpolation to evaluate radial function at distance r.

		Encapsulation required due to autovivification of SmartDict.

		Args:
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoffs
			distance (float): Radial distance from ion
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'

		Returns:
			float: Radial function evaluated at r
		"""

		return self.get_radial(l, zeta).get_value(distance, interpolation=interpolation)

	def get_max_cutoff(self):
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
		position (Vector): 3D Cartesian real space vector for atom position
		bands (SmartDict): Nested dict of complex basis coefficients;
							Accessed by bands[bandEnergy][l][zeta][m]
	"""

	def __init__(self, ionName, atomPos, radials=SmartDict()):
		"""Constructor for atom.

		Args:
			ionName (string): Name of ion (usually from name of .ion file)
			radials (SmartDict, optional): Radial objects, indexed by radialDict[zeta][n][l],
											where all indices are int. Default is empty, radials
											may be added after instantiation
			atomPos (Vector): 3D Cartesian real space vector for atom position
		"""
		Ion.__init__(self, ionName, radials)
		self.position = atomPos
		self.bands = SmartDict()

	def set_ion(self, I):
		"""Copy all attributes from an Ion to this Atom.

		Args:
			I (Ion): Ion object from which to copy attributes
		"""
		self.ionName = I.ionName
		self.radials = I.radials
		self.sortPAOs()

	def within_cutoff_relative(self, relative_position, l=None, zeta=None):
		"""Return true if within cutoff region.

		Args:
			position (Vector): 3D Cartesian real space vector
			l (int, opt.): Orbital angular momentum quantum number; to check specific radial, needs zeta
			zeta (int, opt.): Zeta index; to check specific radial, needs l

		Returns:
			bool: True if within cutoff radius
		"""
		output = False
		distance = abs(relative_position)
		if not (l and zeta):
			if distance <= self.get_max_cutoff():
				output = True
		elif l and zeta:
			if self.has_radial(l, zeta):
				if distance <= self.get_radial(l, zeta).cutoff:
					output = True
		else:
			raise ValueError("Need both l and zeta input, or neither")
		return output

	def within_cutoff(self, position, l=None, zeta=None):
		"""Return true if within cutoff region.

		Args:
			position (Vector): 3D Cartesian real space vector
			l (int, opt.): Orbital angular momentum quantum number; to check specific radial, needs zeta
			zeta (int, opt.): Zeta index; to check specific radial, needs l

		Returns:
			bool: True if within cutoff radius
		"""
		output = False
		relative_position = abs(position - self.position)
		return self.within_cutoff_relative(relative_position, l, zeta)

	def has_coefficient(self, K, E, l, zeta, m):
		"""Check if atom stores coefficient for given orbital.

		Encapsulation required due to autovivification of SmartDict.

		Args:
			K (Vector): 3D Cartesian k-space vector
			E (float): Band energy
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoff but same n and l
			m (int): Azimuthal orbital angular momentum quantum number

		Returns:
			boolean: True if coefficient is stored, false if not
		"""

		output = False
		if K in self.bands:
			if E in self.bands[K]:
				if l in self.bands[K][E]:
					if zeta in self.bands[K][E][l]:
						if m in self.bands[K][E][l][zeta]:
							output = True
		return output

	def add_coefficient(self, K, E, PAO, coefficient):
		"""Add a complex coefficient to self.bands.

		Args:
			K (Vector): 3D Cartesian k-space vector
			E (float): Band energy
			PAO (int): index of PAO as given in .dat file
			coefficient (complex): Coefficient of PAO
		"""

		PAOdata = self.sortPAOs()[PAO - 1]
		l = PAOdata[0]
		zeta = PAOdata[1]
		m = PAOdata[2]

		self.bands[K][E][l][zeta][m] = coefficient

	def get_coefficient(self, K, E, l, zeta, m):
		"""Return complex coefficient for given orbital.

		Args:
			K (Vector): 3D Cartesian k-space vector
			E (float): Band energy
			l (int): Orbital angular momentum quantum number
			zeta (int): Indexes functions with different cutoff but same n and l
			m (int): Azimuthal orbital angular momentum quantum number

		Returns:
			complex: Coefficient for given orbital
		"""
		output = 0.0
		if self.bands.locked:
			output = self.bands[K][E][l][zeta][m]
		elif self.has_coefficient(K, E, l, zeta, m):
			output = self.bands[K][E][l][zeta][m]
		return output

	def get_radial_value_relative(self, l, zeta, relative_position, interpolation='cubic'):
		"""Evaluate radial part of wavefunction

		Args:
			l (int): Orbital angular momentum quantum number
			zeta (int): Zeta index
			position (Vector): 3D Cartesian real space vector
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
		"""
		R = 0.0
		if self.has_radial(l, zeta):
			distance = abs(relative_position)
			R = self.get_radial_value(l, zeta, distance, interpolation=interpolation)
		return R

	def get_sph_relative(self, l, m, position):
		"""Evaluate spherical harmonic with atom as origin

		Args:
			l (int): Orbital angular momentum quantum number
			m (int): Azimuthal orbital angular momentum quantum number
			position (Vector): 3D Cartesian real space vector
		"""
		position2 = position - self.position
		return sph(l, m, position2)

	def get_basis_point(self, l, zeta, m, position, interpolation='cubic'):
		"""Evaluate basis function (radial part * spherical harmonic)"""
		R = self.get_radial_value_relative(l, zeta, position - self.position, interpolation=interpolation)
		Y = 0.0
		if R != 0:  # If R == 0, basis point is 0, no need to calculate Y
			Y = self.get_sph_relative(l, m, position - self.position)
		return R*Y

	def get_psi(self, K, E, position, interpolation='cubic', basisPoint=SmartDict(), local=False):
		"""Evaluate wavefunction contribution from this atom.

		Args:
			K (Vector): 3D Cartesian k-space vector
			E (float): Band energy
			position (Vector): 3D Cartesian real space vector
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
			basisPoint (SmartDict, opt.): Basis function values indexed by [l][zeta][m]
			local (bool, opt.): If true, locality has already been checked ie. point is known to lie within atom cutoff

		Returns:
			complex: Wavefunction value
		"""
		psi = complex(0.0, 0.0)

		if not local:
			if not self.within_cutoff(position):
				return psi

		if not basisPoint:
			for l in self.radials:
				for zeta in self.radials[l]:
					for m in range(l, l+1):
						basisPoint[l][zeta][m] = self.get_basis_point(l, zeta, m, position, interpolation=interpolation)
		for l in basisPoint:
			for zeta in basisPoint[l]:
				for m in basisPoint[l][zeta]:
					if basisPoint[l][zeta][m] != 0:
						coefficient = self.get_coefficient(K, E, l, zeta, m)
						psi += basisPoint[l][zeta][m] * coefficient
		return psi

	def apply_factor(self, factor, K, E):
		"""Apply a normalisation factor to all coefficients for a given energy and k-point.

		Args:
			factor (float): Normalisation factor to be applied
			K (Vector): 3D Cartesian k-space vector
			E (float): Energy of band to which to apply factor
		"""
		# Loop over all coefficients
		for l in self.bands[K][E]:
			for zeta in self.bands[K][E][l]:
				for m in self.bands[K][E][l][zeta]:
					# Apply factor
					self.bands[K][E][l][zeta][m] *= factor
