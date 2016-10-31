
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from vector import Vector

class Radial:

	"""Stores the radial part of basis function and metadata,
	ie. quantum numbers (n and l) and zeta index."""

	def __init__(self, zeta, n, l, r, R):
		self.zeta = zeta
		self.n = n
		self.l = l
		self.r = r
		self.R = R
	
class Ion(object):

	"""Stores Radial objects."""

	# Normalisation factors for spherical harmonic
	# Numbers in variable names refer to l and abs(m)
	sph00 = 1.0/2.0 * np.sqrt(1.0   / math.pi)

	sph11 = 1.0/2.0 * np.sqrt(3.0   / (2.0*math.pi))
	sph10 = 1.0/2.0 * np.sqrt(3.0   / math.pi)

	sph22 = 1.0/4.0 * np.sqrt(15.0  / (2.0*math.pi)) 
	sph21 = 1.0/2.0 * np.sqrt(15.0  / (2.0*math.pi))
	sph20 = 1.0/4.0 * np.sqrt(5.0   / math.pi)

	sph33 = 1.0/8.0 * np.sqrt(35.0  / math.pi)
	sph32 = 1.0/4.0 * np.sqrt(105.0 / (2.0*math.pi))
	sph31 = 1.0/8.0 * np.sqrt(21.0  / math.pi)
	sph30 = 1.0/4.0 * np.sqrt(7.0   / math.pi)

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
		"""Find closest r value in Radial to given r and return R(r)"""
		Rad = self.Rads[zeta][n][l]
		rvalues = Rad.r
		Rvalues = Rad.R

		bestIndex = 0
		lowestdif = 1000000
		for i in range(0, len(rvalues)):
			dif = r - rvalues[i]
			if abs(lowestdif) > abs(dif):
				lowestdif = dif
				bestIndex = i
		if lowestdif > 0.5:
			return 0.0
		else:
			return Rvalues[bestIndex]

	@classmethod
	def sph(cls, l, m, x, y, z):
		"""Return cartesian spherical harmonic.
		Valid indices are l > 0; -l < m < +l.
		Currently, only l <= 3 is supported.
		
		Input:
		l:			Orbital angular momentum quantum number
		m:			Azimuthal quantum number
		x, y, z:	Cartesian coordinates

		Output:
		harm:		Spherical harmonic for given l and m, evaluated at (x, y, z)"""

		# Normalise coordinates to unit magnitude
		# If at origin, return 0
		magnitude = np.sqrt(x**2 + y**2 + z**2)
		if magnitude != 0:
			x = x / magnitude
			y = y / magnitude
			z = z / magnitude
		else:
			return 0.0

		# Choose correct form of spherical harmonic depending on l and m
		harm = 0.0
		if l == 0:
			harm =      cls.sph00
		elif l == 1:
			if m == 1:
				harm =  cls.sph11 * x
			elif m == -1:
				harm =  cls.sph11 * y
			elif m == 0:
				harm =  cls.sph10 * z
		elif l == 2:
			if m == 2:
				harm =  cls.sph22 * (x**2 - y**2)
			elif m == -2:
				harm =  cls.sph22 * (x**2 - y**2)
			elif m == 1:
				harm =  cls.sph21 * x*z
			elif m == -1:
				harm = -cls.sph21 * x*z
			elif m == 0:
				harm =  cls.sph20 * (3*z**2 - 1)
		elif l == 3:
			if m == 3:
				harm =  -cls.sph33 * (x**3 - 3*x*y**2)
			elif m == -3:
				harm =   cls.sph33 * (x**3 - 3*x*y**2)
			elif m == 2:
				harm = 	 cls.sph32 * (x**2 - y**2)*z
			elif m == -2:
				harm =   cls.sph32 * (x**2 - y**2)*z
			elif m == 1:
				harm =  -cls.sph31 * (5*z**2 - 1)*x
			elif m == -1:
				harm =   cls.sph31 * (5*z**2 - 1)*x
			elif m == 0:
				harm =   cls.sph30 * (5*z**2 - 3)*z
		return harm

	def plotBasis(self, zeta, n, l, m, axis, minimum=-8, maximum=8, planeValue=0.00001, step=0.1):
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
		space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]
		Y = np.empty_like(space1)
		R = np.empty_like(space1)
		psi = np.empty_like(space1)

		maxpsi = 0.1
		with PdfPages('pdfs/'+plotname+'.pdf') as pdf:
			print 'Creating '+plotname+'.pdf'
			for i in range(0, int((maximum-minimum)/step)):
				for j in range(0, int((maximum-minimum)/step)):
					if axis == 'z':
						Y[i,j] = self.sph(l,m,space2[i,j],space1[i,j],planeValue)
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$y$ / $a_0$')
					if axis == 'y':
						Y[i,j] = self.sph(l,m,space2[i,j],planeValue,space1[i,j])
						plt.xlabel('$x$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')
					if axis == 'x':
						Y[i,j] = self.sph(l,m,planeValue,space2[i,j],space1[i,j])
						plt.xlabel('$y$ / $a_0$')
						plt.ylabel('$z$ / $a_0$')
					
					R[i, j] = self.getRadialValue(zeta, n, l, np.sqrt(space1[i,j]**2+space2[i,j]**2+planeValue**2))
					psi[i,j] = Y[i,j]*R[i,j]
					#if space1[i, j] == 5.0 and space2[i, j] == 2.0:
					#	psi[i,j] = 5
					if abs(psi[i,j]) > maxpsi:
						maxpsi = abs(psi[i,j])

			#plt.contour(space2, space1, psi, 8, cmap=plt.cm.seismic)
			plt.imshow(psi, interpolation='bilinear', origin='center', cmap=plt.cm.bwr, extent=(minimum,maximum,minimum,maximum),vmin=-maxpsi, vmax=maxpsi)
			plt.colorbar()
			plt.grid()
			axes = ['x', 'y', 'z']
			axes.remove(axis)
			ttl = self.name+' Basis Function for \n \n $\zeta='+str(zeta)+'$, $n='+str(n)+'$, $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+axes[0]+'-'+axes[1]+'$ plane'
			plt.title(ttl)
			pdf.savefig()
			plt.close()

	@classmethod
	def plotSPH(cls, l, m):
		minimum = -8
		maximum = 8

		step = 0.1
		x,y = np.mgrid[minimum:maximum:step, minimum:maximum:step]
		
		Y = np.empty_like(x)
		maxY = 0
		for i in range(0, int((maximum-minimum)/step)):
			for j in range(0, int((maximum-minimum)/step)):
				Y[i, j] = cls.sph(l,m,x[i,j],y[i,j],1)
				if abs(Y[i,j]) > maxY:
					maxY = abs(Y[i,j])

		plt.contour(x, y, Y, 0)
		plt.imshow(Y, interpolation='bilinear', cmap=plt.cm.seismic, extent=(minimum,maximum,minimum,maximum),vmin=-maxY, vmax=maxY)
		plt.colorbar()
		ttl = 'Spherical Harmonic for $l='+str(l)+'$, $m_l='+str(m)+'$'
		plt.title(ttl)
		plt.show()


