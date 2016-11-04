
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from packages.sph import sph
from packages.smartDict import SmartDict

class Radial(object):

	"""Stores the radial part of basis function and metadata,
	ie. quantum numbers (n and l) and zeta index.
	"""

	def __init__(self, zeta, n, l, r, R, cutoff):
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
		self.r = r
		self.R = R
		self.cutoff = cutoff

	def getValue(self, distance):
		"""Use linear interpolation to evaluate radial function at distance.
		
		Args:
		    distance (double): Distance from origin in Bohr radii
		
		Returns:
		    double: Value of radial part
		"""
		if distance > self.cutoff or distance > self.r[-1]:
			value = 0.0
		else:
			# Find the first r value larger than distance
			i = 0
			while self.r[i] < distance:
				i = i + 1

			# Get nearest stored values
			x1 = self.r[i-1]
			x2 = self.r[i]
			y1 = self.R[i-1]
			y2 = self.R[i]

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

		# Initialise dict entry
		# if not self.radials.has_key(zeta):
		# 	self.radials[zeta] = {}
		# if not self.radials[zeta].has_key(n):
		# 	self.radials[zeta][n] = {}
		# if not l in self.nl[n]:
		# 	self.nl[n].append(l)

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

	def plotBasis(self, zeta, n, l, m, axis, minimum=-8, maximum=8, planeValue=0.0, step=0.1):
		"""Plots cross-section of basis function of ion to pdf.
		All lengths measured in bohr radii (a0).

		Input:
		zeta:		Zeta index of Radial
		n:			Principal quantum number for Radial
		l:			Orbital angular momentum quantum number for Radial and spherical harmonic
		m:			Azimuthal quantum number for spherical harmonic
		axis:		Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum:	Minimum value of coordinates measured in a0; Default is -8
		maximum:	Maximum value of coordinates measured in a0; Default is +8
		planeValue:	Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
		step:		Interval between Cartesian mgrid points, measured in a0; Default is 0.1"""

		plotname = 'Basis_'+self.ionName+'_'+str(zeta)+'_'+str(n)+'_'+str(l)+'_'+str(m)+'_'+axis

		# Initialise meshes
		# 2D cartesian mesh (x, y, or z axis determined later)
		space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]

		Y = np.empty_like(space1) # Spherical Harmonic mesh
		R = np.empty_like(space1) # Radial Function mesh
		psi = np.empty_like(space1) # Basis Function mesh (psi = R*Y)

		maxpsi = 0.1 # Colour plot sets limits to -maxpsi to +maxpsi

		# Plot functions to pdf
		with PdfPages('pdfs/' + plotname + '.pdf') as pdf:
			# Loop over all mesh points
			for i in range(0, int((maximum - minimum) / step)):
				for j in range(0, int((maximum - minimum) / step)):
					# Use axis variable to determine which axes space1 and space2 refer to
					# Evaluate spherical harmonic at mesh point
					if axis == 'z':
						Y[i, j] = sph(l, m, space2[i, j], space1[i, j], planeValue)
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$y$ / $a_0$')
					if axis == 'y':
						Y[i, j] = sph(l, m, space2[i, j], planeValue, space1[i, j])
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')
					if axis == 'x':
						Y[i, j] = sph(l, m, planeValue, space2[i, j], space1[i, j])
						plt.xlabel('$y$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')

					# Estimate value of Radial for mesh point and get psi
					R[i, j] = self.getRadialValue(zeta, n, l, np.sqrt(space1[i, j]**2 +
						                                                 space2[i, j]**2 +
						                                                 planeValue**2))
					psi[i, j] = Y[i, j] * R[i, j]

					# Update maxpsi
					if abs(psi[i, j]) > maxpsi:
						maxpsi = abs(psi[i, j])

			# Setup plot
			plt.imshow(psi, interpolation='bilinear', origin='center', cmap=plt.cm.bwr,
				          extent=(minimum, maximum, minimum, maximum), vmin=-maxpsi, vmax=maxpsi)
			plt.colorbar()
			plt.grid()
			axes = ['x', 'y', 'z']
			axes.remove(axis)
			ttl = (self.ionName+' Basis Function for \n \n $\zeta='+str(zeta)+'$, $n='+
				      str(n)+'$, $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+axes[0]+'-'+axes[1]+'$ plane')
			plt.title(ttl)

			# Save to pdf
			pdf.savefig()
			plt.close()
			print 'Finished '+plotname+'.pdf'

	@classmethod
	def plotSPH(cls, l, m, axis, minimum=-8, maximum=8, planeValue=0.0, step=0.1):
		"""Plots cross-section of spherical harmonic to pdf.
		All lengths measured in bohr radii (a0).

		Input:
		l:			Orbital angular momentum quantum number for spherical harmonic
		m:			Azimuthal quantum number for spherical harmonic
		axis:		Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum:	Minimum value of coordinates measured in a0; Default is -8
		maximum:	Maximum value of coordinates measured in a0; Default is +8
		planeValue:	Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
		step:		Interval between Cartesian mgrid points, measured in a0; Default is 0.1"""

		plotname = 'SPH_'+str(l)+'_'+str(m)+'_'+axis

		# Initialise meshes
		# 2D cartesian mesh (x, y, or z axis determined later)
		space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]
		Y = np.empty_like(space1) # Spherical Harmonic mesh

		maxY = 0.1 # Colour plot sets limits to -maxY and +maxY
		with PdfPages('pdfs/'+plotname+'.pdf') as pdf:
			for i in range(0, int((maximum - minimum) / step)):
				for j in range(0, int((maximum - minimum) / step)):
					# Use axis variable to determine which axes space1 and space2 refer to
					# Evaluate spherical harmonic at mesh point
					if axis == 'z':
						Y[i, j] = sph(l, m, space2[i, j], space1[i, j], planeValue)
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$y$ / $a_0$')
					if axis == 'y':
						Y[i, j] = sph(l, m, space2[i, j], planeValue, space1[i, j])
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')
					if axis == 'x':
						Y[i, j] = sph(l, m, planeValue, space2[i, j], space1[i, j])
						plt.xlabel('$y$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')

					# Update maxY
					if abs(Y[i, j]) > maxY:
						maxY = abs(Y[i, j])

			# Setup plot
			plt.imshow(Y, interpolation='bilinear', origin='center', cmap=plt.cm.bwr,
				          extent=(minimum, maximum, minimum, maximum), vmin=-maxY, vmax=maxY)
			plt.colorbar()
			plt.grid()
			axes = ['x', 'y', 'z']
			axes.remove(axis)
			ttl = ('Spherical Harmonic for \n \n $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+
				      axes[0]+'-'+axes[1]+'$ plane')
			plt.title(ttl)

			# Save to pdf
			pdf.savefig()
			plt.close()
			print 'Finished '+plotname+'.pdf'

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
		PAOdata = self.sortPAOs()[PAO-1]
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

		# If r is beyond atoms range, return 0.0
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
