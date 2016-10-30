
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from vector import Vector

class Radial:
	"""Stores the radial part of basis function and metadata, ie.
	quantum numbers (n and l) and zeta index"""

	def __init__(self, zeta, n, l, r, R):
		self.zeta = zeta
		self.n = n
		self.l = l
		self.r = r
		self.R = R
	
class Ion(object):
	"""Stores Radial objects"""
	sph00 = 1.0 / (2.0 * np.sqrt(math.pi))
	sph11 = 1.0/2.0 * np.sqrt(3.0 / (2.0*math.pi))
	sph10 = 1.0/2.0 * np.sqrt(3.0 / math.pi)
	sph22 = 1.0/4.0 * np.sqrt(15.0 / (2.0*math.pi)) 
	sph21 = 1.0/2.0 * np.sqrt(15.0 / (2.0*math.pi))
	sph20 = 1.0/4.0 * np.sqrt(5.0 / math.pi)

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

		harm = 0.0
		if l == 0:
			harm = cls.sph00
		elif l == 1:
			if m == 1:
				harm = cls.sph11 * x
			elif m == -1:
				harm = cls.sph11 * y
			elif m == 0:
				harm =  cls.sph10 * z
		elif l == 2:
			if m == 2:
				harm =  cls.sph22 * (x**2 - y**2)
			elif m == -2:
				harm =  cls.sph22 * (x**2 - y**2)
			elif m == 1:
				harm = -cls.sph21 * x*z
			elif m == -1:
				harm =  cls.sph21 * x*z
			elif m == 0:
				harm =  cls.sph20 * 3*z**2
		return harm

	def plotBasis(self, zeta, n, l, m, axis):
		plotname = 'SPH_'+str(l)+str(m)

		minimum = -8
		maximum = 8

		step = 0.1
		space = np.mgrid[minimum:maximum:step, minimum:maximum:step]
		Y = self.sumInQuad(space[0],space[1])
		R = self.sumInQuad(space[0],space[1])
		psi = self.sumInQuad(space[0],space[1])
		maxpsi = 0
		for i in range(0, int((maximum-minimum)/step)):
			for j in range(0, int((maximum-minimum)/step)):
				if axis == 'z':
					Y[i,j] = self.sph(l,m,space[0][i,j],space[1][i,j],0.01)
					plt.xlabel('$x$ / $a_0$')
					plt.ylabel('$y$ / $a_0$')
				if axis == 'y':
					Y[i,j] = self.sph(l,m,space[0][i,j],0.01,space[1][i,j])
					plt.xlabel('$x$ / $a_0$')
					plt.ylabel('$z$ / $a_0$')
				if axis == 'x':
					Y[i,j] = self.sph(l,m,0.01,space[0][i,j],space[1][i,j])
					plt.xlabel('$y$ / $a_0$')
					plt.ylabel('$z$ / $a_0$')
				
				R[i, j] = self.getRadialValue(zeta, n, l, np.sqrt(space[0][i,j]**2+space[1][i,j]**2+1))
				psi[i,j] = Y[i,j]*R[i,j]
				if abs(psi[i,j]) > maxpsi:
					maxpsi = abs(psi[i,j])


		plt.imshow(psi, interpolation='bilinear', cmap=plt.cm.seismic, extent=(minimum,maximum,minimum,maximum),vmin=-maxpsi, vmax=maxpsi)
		plt.colorbar()
		axes = ['x', 'y', 'z']
		axes.remove(axis)
		ttl = 'Basis Function, $\zeta='+str(zeta)+'$, $n='+str(n)+'$, $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+axes[0]+'-'+axes[1]+'$ plane'
		plt.title(ttl)
		plt.show()

	@staticmethod
	def sumInQuad(x, y):
		return x**2 + y**2

	@classmethod
	def plotSPH(cls, l, m):
		minimum = -8
		maximum = 8

		step = 0.1
		x,y = np.mgrid[minimum:maximum:step, minimum:maximum:step]
		
		Y = cls.sumInQuad(x,y)
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


