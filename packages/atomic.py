
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
		"""Constructs radial part of a basis function
		
		Args:
		    zeta (int): Indexes functions with the same n and l, but different cutoff
		    n (int): Principal quantum number
		    l (int): Orbital angular momentum quantum number
		    r (float[]): List of radial distance values; Measured in Bohr radii
		    R (float[]): List of radial function values; Same length as r
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
		    distance (double): Distance from origin in Bohr radii
		
		Returns:
		    double: Value of radial part
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
	Stores data on an ion represented using a certain basis.

	Attributes:
	ionName (string): Name of ion (usually from name of .ion file)
	radials (nested dict): Radial objects stored in nested dict;
							Indexed by radials[zeta][n][l]
	orbitals (nested dict): Which 

	"""
	def __init__(self, ionName, radialDict=SmartDict()):
		"""Summary
		
		Args:
		    ionName (string): Name of ion (usually from name of .ion file)
		    radialDict (int, dict()): Description
		"""
		self.ionName = ionName
		self.radials = radialDict # Radial objects; accessed by self.radials[zeta][n][l]

		# # Loop over all zeta indices in radialDict
		# for zeta, nDict in radialDict:
		# 	# If zeta not already added to orbitalDict, create empty dict
		# 	if not orbitalDict.has_key(zeta):
		# 		orbitalDict[zeta] = {}
		# 	# Loop over all n for this zeta
		# 	for n, lDict in nDict:
		# 		# If n not already added to orbitalDict, create empty dict
		# 		if not orbitalDict[zeta].has_key(n):
		# 			orbitalDict[zeta][n] = []
		# 		# Add l values, lowest to highest, to orbitalDict
		# 		for l in sorted(lDict):
		# 			orbitalDict[zeta][n].append(l)
		# self.orbitals = orbitalDict


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
		dict if given an invalid key."""
		output = False
		if self.radials.has_key(zeta):
			if self.radials[zeta].has_key(n):
				if self.radials[zeta][n].has_key(l):
					return True
		return output

	def addRadial(self, radial):
		"""Adds radial to self.radials. Overwrites radial with same metadata (zeta, n, l)"""
		# Get metadata
		zeta = radial.zeta
		n = radial.n
		l = radial.l

		# Add Radial
		self.radials[zeta][n][l] = radial
		self.sortPAOs()

	def getRadial(self, zeta, n, l):
		output = None
		if self.hasRadial(zeta, n, l):
			output = self.radials[zeta][n][l]
		return output

	def getRadialValue(self, zeta, n, l, r):
		"""Use linear interpolation to evaluate R at r."""
		output = None
		if self.hasRadial(zeta, n, l):
			output = self.radials[zeta][n][l].getValue(r)
		return output

	def getMaxCutoff(self):
		"""Return the maximum cutoff of radius of Radial bases, beyond which
		the radial part is defined to be 0."""
		maxcut = 0.0
		for zeta in self.radials:
			for n in self.radials[zeta]:
				for l in self.radials[zeta][n]:
					if maxcut < self.radials[zeta][n][l].cutoff:
						maxcut = self.radials[zeta][n][l].cutoff
		return maxcut

class Atom(Ion):

	"""Stores information about an atom, primarily the Ions and coefficients.
	UNFINISHED"""

	def __init__(self, name, x, y, z):
		Ion.__init__(self, name)
		self.x = x
		self.y = y
		self.z = z
		self.coeffs = SmartDict()

	def setIon(self, I):
		'''Copy all attributes from an Ion to this Atom'''
		self.ionName = I.ionName
		self.radials = I.radials
		self.sortPAOs()

	def hasCoeff(self, zeta, n, l, m):

		output = False
		if self.coeffs.has_key(zeta):
			if self.coeffs[zeta].has_key(n):
				if self.coeffs[zeta][n].has_key(l):
					if self.coeffs[zeta][n][l].has_key(m):
						output = True
		return output

	def addCoeff(self, PAO, coeff):
		'''Add a complex coefficient to self.coeffs'''

		# Get zeta, n, l, m from the PAO index as given in .dat file
		PAOdata = self.sortPAOs()[PAO - 1]
		zeta = PAOdata[0]
		n = PAOdata[1]
		l = PAOdata[2]
		m = PAOdata[3]

		self.coeffs[zeta][n][l][m] = coeff

	def getCoeff(self, zeta, n, l, m):
		output = None
		if self.hasCoeff(zeta, n, l, m):
			output = self.coeffs[zeta][n][l][m]
		return output

	def getPsi(self, x, y, z):
		"""Evaluate the wavefuncion at (x,y,z) due to this atom only"""

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
							coeff = self.coeffs[zeta][n][l][m]

							# Calculate and add contribution of basis function
							psiReal = R*Y*coeff.real
							psiImag = R*Y*coeff.imag
							psi += complex(psiReal, psiImag)
		return psi
