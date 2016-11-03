import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from packages.sph import sph

class Radial(object):

	"""Stores the radial part of basis function and metadata,
	ie. quantum numbers (n and l) and zeta index."""

	def __init__(self, zeta, n, l, r, R, cutoff):
		self.zeta = zeta
		self.n = n
		self.l = l
		self.r = r
		self.R = R
		self.cutoff = cutoff

	def getValue(self, distance):
		"""Use linear interpolation to evaluate the radial function
		at distance."""

		if r > self.cutoff:
			return 0.0
		else:
			i = 0
			r0 = 0.0
			r1 = 0.0
			while self.r[i] < distance:
				i = i + 1

			r0 = self.r[i-1]
			r1 = self.r[i]
			R0 = self.r[i-1]
			R1 = self.r[i]

			value = R0 + (r - r0) * (R1 - R0) / (r1 - r0)
			return value

class Ion(object):

	"""Stores information about an ion, primarily its basis functions."""

	def __init__(self, name, radialDict=None):
		self.ionName = name

		#self.zetas = 1 # 1 = SZ; 2 = DZ; 3 = TZ
		self.zetas = 1
		self.nl = {}
		if radialDict:
			self.radials = radialDict # Radial objects; accessed by self.radials[zeta][n][l]
			for zeta in self.radials.keys():
				if zeta > self.zetas:
					self.zetas = zeta
				for n in self.radials[zeta]:
					if n not in self.nl.keys():
						self.radials[zeta][n] = []
					for l in self.radials[zeta][n]:
						if l not in self.nl[n]:
							self.nl[n].append(l)
		else:
			self.radials = {}

	def sortPAOs(self):
		sortedPAOs = []
		for zeta in range(1, self.zetas+1):
			nList = sorted(self.nl.keys())
			for n in nList:
				lList = sorted(self.nl[n])
				for l in lList:
					for m in range(-l, l+1):
						sortedPAOs.append([zeta, n, l, m])
		return sortedPAOs

	def addRadial(self, radial):
		"""Adds Radial to self.radials. Overwrites radial with same metadata (zeta, n, l)"""
		# Get metadata
		zeta = radial.zeta
		n = radial.n
		l = radial.l

		# Set max value of zeta
		if zeta > self.zetas:
			self.zetas = zeta

		# Initialise dict entry
		if not self.radials.has_key(zeta):
			self.radials[zeta] = {}
		if not self.radials[zeta].has_key(n):
			self.radials[zeta][n] = {}
			self.nl[n] = []
		if not l in self.nl[n]:
			self.nl[n].append(l)

		# Add Radial
		self.radials[zeta][n][l] = radial
		self.sortPAOs()

	def getRadial(self, zeta, n, l):
		return self.radials[zeta][n][l]

	def getRadialValue(self, zeta, n, l, r):
		"""Use linear interpolation to evaluate R at r."""
		return self.radials[zeta][n][l].getValue(r)

	def getMaxCutoff(self):
		"""Return the maximum cutoff of radius of Radial bases, beyond which
		the radial part is defined to be 0."""
		maxcut = 0.0
		for zeta in range(1, self.zetas+1):
			for n in self.radials[zeta].keys():
				for l in self.radials[zeta][n].keys():
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
		self.coeffs = {}

	def setIon(self, I):
		'''Copy all attributes from an Ion to this Atom'''
		self.radials = I.radials
		self.zetas = I.zetas
		self.nl = I.nl
		self.sortPAOs()

	def addCoeff(self, PAO, coeff):
		'''Add a complex coefficient to self.coeffs'''

		# Get zeta, n, l, m from the PAO index as given in .dat file
		PAOdata = self.sortPAOs()[PAO-1]
		zeta = PAOdata[0]
		n = PAOdata[1]
		l = PAOdata[2]
		m = PAOdata[3]

		# Initialise dict entry
		if not self.coeffs.has_key(zeta):
			self.coeffs[zeta] = {}
		if not self.coeffs[zeta].has_key(n):
			self.coeffs[zeta][n] = {}
		if not l in self.coeffs[zeta][n]:
			self.coeffs[zeta][n][l] = {}

		# Set the coefficient
		self.coeffs[zeta][n][l][m] = coeff

	def getPsi(self, x, y, z):
		"""Evaluate the wavefuncion at (x,y,z)"""

		# Get coords relative to atom
		relx = x - self.x
		rely = y - self.y
		relz = z - self.z
		# Get relative distance to atom
		r = np.sqrt(relx**2 + rely**2 + relz**2)

		if r > self.getMaxCutoff():
			return complex(0.0, 0.0)
		else:
			psi = complex(0.0, 0.0)
			for zeta in range(1, self.zetas+1):
				for n in self.nl.keys():
					for l in self.nl[n]:
						R = self.getRadialValue(zeta, n, l, r)
						for m in range(-l, l+1):
							Y = sph(l, m, relx, rely, relz)
							coeff = self.coeffs[zeta][n][l][m]
							psiReal = R*Y*coeff.real
							psiImag = R*Y*coeff.imag
							psi += complex(psiReal, psiImag)
							#print psi
			return psi
