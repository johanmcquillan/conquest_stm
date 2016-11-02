
from sph import *
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class Radial:

	"""Stores the radial part of basis function and metadata,
	ie. quantum numbers (n and l) and zeta index."""

	def __init__(self, zeta, n, l, r, R, cutoff):
		self.zeta = zeta
		self.n = n
		self.l = l
		self.r = r
		self.R = R
		self.cutoff = cutoff
	
class Ion(object):

	"""Stores information about an ion, primarily its basis functions."""

	def __init__(self, name):
		self.name = name
		self.nl = {}
		self.zetas = 1 # 1 = SZ; 2 = DZ; 3 = TZ
		self.Rads = {} # Radial objects; accessed by self.Rads[zeta][n][l]

	def addRadial(self, radial):
		"""Adds Radial to self.Rads. Overwrites radial with same metadata"""
		# Get metadata
		zeta = radial.zeta
		n = radial.n
		l = radial.l

		# Set max value of zeta
		if zeta > self.zetas:
			self.zetas = zeta

		# Initialise dict entry
		if not self.Rads.has_key(zeta):
			self.Rads[zeta] = {}
		if not self.Rads[zeta].has_key(n):
			self.Rads[zeta][n] = {}
			self.nl[n] = []
		if not l in self.nl[n]:
			self.nl[n].append(l)

		# Add Radial
		self.Rads[zeta][n][l] = radial

	def getRadial(self, zeta, n, l):
		return self.Rads[zeta][n][l]

	def getRadialValue(self, zeta, n, l, r):
		"""Use linear interpolation to evaluate R at r."""

		# Get data
		Rad = self.Rads[zeta][n][l]
		if r > Rad.cutoff:
			return 0.0
		else:
			rvalues = Rad.r
			Rvalues = Rad.R

			i = 0
			r0 = 0.0
			r1 = 0.0
			while rvalues[i] < r:
				i = i + 1

			r0 = rvalues[i-1]
			r1 = rvalues[i]
			R0 = Rvalues[i-1]
			R1 = Rvalues[i]

			R = R0 + (r - r0) * (R1 - R0) / (r1 - r0)
			return R

	def getMaxCutoff(self):
		"""Return the maximum cutoff of radius of Radial bases, beyond which 
		the radial part is defined to be 0."""
		maxcut = 0.0
		for z in range(1, self.zetas+1):
			for n in self.Rads[zeta].keys():
				for l in self.Rads[zeta][n].keys():
					if maxcut < self.Rads[zeta][n][l].cutoff:
						maxcut = self.Rads[zeta][n][l].cutoff
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

		plotname = 'Basis_'+self.name+'_'+str(zeta)+'_'+str(n)+'_'+str(l)+'_'+str(m)+'_'+axis

		# Initialise meshes
		space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step] # 2D cartesian mesh (x, y, or z axis determined later)
		Y = np.empty_like(space1) # Spherical Harmonic mesh
		R = np.empty_like(space1) # Radial Function mesh
		psi = np.empty_like(space1) # Basis Function mesh (psi = R*Y)

		maxpsi = 0.1 # Colour plot sets limits to -maxpsi to +maxpsi

		# Plot functions to pdf
		with PdfPages('pdfs/'+plotname+'.pdf') as pdf:
			# Loop over all mesh points
			for i in range(0, int((maximum-minimum)/step)):
				for j in range(0, int((maximum-minimum)/step)):
					# Use axis variable to determine which axes space1 and space2 refer to
					# Evaluate spherical harmonic at mesh point
					if axis == 'z':
						Y[i,j] = sph(l,m,space2[i,j],space1[i,j],planeValue)
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$y$ / $a_0$')
					if axis == 'y':
						Y[i,j] = sph(l,m,space2[i,j],planeValue,space1[i,j])
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')
					if axis == 'x':
						Y[i,j] = sph(l,m,planeValue,space2[i,j],space1[i,j])
						plt.xlabel('$y$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')
					
					# Estimate value of Radial for mesh point and get psi
					R[i, j] = self.getRadialValue(zeta, n, l, np.sqrt(space1[i,j]**2+space2[i,j]**2+planeValue**2))
					psi[i,j] = Y[i,j]*R[i,j]

					# Update maxpsi
					if abs(psi[i,j]) > maxpsi:
						maxpsi = abs(psi[i,j])

			# Setup plot
			plt.imshow(psi, interpolation='bilinear', origin='center', cmap=plt.cm.bwr, extent=(minimum,maximum,minimum,maximum),vmin=-maxpsi, vmax=maxpsi)
			plt.colorbar()
			plt.grid()
			axes = ['x', 'y', 'z']
			axes.remove(axis)
			ttl = self.name+' Basis Function for \n \n $\zeta='+str(zeta)+'$, $n='+str(n)+'$, $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+axes[0]+'-'+axes[1]+'$ plane'
			plt.title(ttl)

			# Save to pdf
			pdf.savefig()
			plt.close()
			print 'Finished '+plotname+'.pdf'

	@classmethod
	def plotSPH(cls, l, m, axis, minimum=-8, maximum=8, planeValue=0.00001, step=0.1):
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
		space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step] # 2D cartesian mesh (x, y, or z axis determined later)
		Y = np.empty_like(space1) # Spherical Harmonic mesh

		maxY = 0.1 # Colour plot sets limits to -maxY and +maxY
		with PdfPages('pdfs/'+plotname+'.pdf') as pdf:
			print 'Creating '+plotname+'.pdf'
			for i in range(0, int((maximum-minimum)/step)):
				for j in range(0, int((maximum-minimum)/step)):
					# Use axis variable to determine which axes space1 and space2 refer to
					# Evaluate spherical harmonic at mesh point
					if axis == 'z':
						Y[i,j] = sph(l,m,space2[i,j],space1[i,j],planeValue)
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$y$ / $a_0$')
					if axis == 'y':
						Y[i,j] = sph(l,m,space2[i,j],planeValue,space1[i,j])
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')
					if axis == 'x':
						Y[i,j] = sph(l,m,planeValue,space2[i,j],space1[i,j])
						plt.xlabel('$y$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')
					
					# Update maxY
					if abs(Y[i,j]) > maxY:
						maxY = abs(Y[i,j])

			# Setup plot
			plt.imshow(Y, interpolation='bilinear', origin='center', cmap=plt.cm.bwr, extent=(minimum,maximum,minimum,maximum),vmin=-maxY, vmax=maxY)
			plt.colorbar()
			plt.grid()
			axes = ['x', 'y', 'z']
			axes.remove(axis)
			ttl = 'Spherical Harmonic for \n \n $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+axes[0]+'-'+axes[1]+'$ plane'
			plt.title(ttl)
			
			# Save to pdf
			pdf.savefig()
			plt.close()
			print 'Finished '+plotname+'.pdf'

class Atom(Ion):

	"""Stores information about an atom, primarily the Ion and basis coefficients.
	UNFINISHED"""

	def __init__(self, name):
		Ion.__init__(self, name)





