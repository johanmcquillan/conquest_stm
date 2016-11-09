
import cmath
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm, colors

from packages.sph import sph
from packages.smartDict import SmartDict

class Plotter(object):

	"""Stores dict of Ion objects and provides methods to plot data to pdf.

	Attributes:
		fname (string): Name prefix for output .pdf files
		ions (dict): Stores Ion objects indexed by ionName
	"""

	# Spectroscopic notation dictionary
	spectral = {0 : 's', 1 : 'p', 2 : 'd', 3 : 'f'}

	def __init__(self, filename, ions):
		self.fname = filename
		self.ions = ions

	def plotRadials(self, points=500, printStatus=False, spectro=True):
		"""Plot all radial functions from self.ions to 'self.filename'_radials.pdf

		Args:
			points (int, optional): Number of points for plot
			printStatus (boolean, optional): If true, print notification when finishing a plot
			spectro (boolean, optional): If true, use spectroscopic notation
		"""
		with PdfPages('pdfs/'+self.fname+'_radials.pdf') as pdf:

			# Plot all functions for the same ion on one graph
			ionNames = sorted(self.ions.keys())
			for ionName in ionNames:
				ion = self.ions[ionName]

				# Setup plot
				plt.title('PAOs for '+ionName+'.ion')
				plt.xlabel('Radial Distance, $r$ / $a_0$')
				plt.ylabel('$R_{nl}(r)$')
				plt.minorticks_on()
				plt.grid(b=True, which='major', alpha=0.45, linestyle='-')
				plt.grid(b=True, which='minor', alpha=0.10, linestyle='-')

				# Loop over all radial functions

				for zeta in ion.radials:
					for n in ion.radials[zeta]:
						for l in ion.radials[zeta][n]:
							# Get Radial and data from ion
							radial = ion.getRadial(zeta, n, l)
							step = radial.cutoff / points
							r = np.arange(0.0, radial.cutoff, step)
							R = np.empty_like(r)
							for i in range(0, len(r)):
								R[i] = radial.getValue(r[i])

							# Add radial info to legend and add to plot
							# If spectro, use spectroscopic notation for legend
							if spectro:
								label = '$\zeta ='+str(zeta)+'$, $'+str(n)+self.spectral[l]+'$'
							else:
								label = '$\zeta ='+str(zeta)+'$, $n='+str(n)+'$, $l='+str(l)+'$'
							plt.plot(r, R, label=label)

							if printStatus:
								print "Finished Radial "+ion.ionName+"_"+str(zeta)+"_"+str(n)+"_"+str(l)

				# Add plot to pdf and reset plt
				plt.legend()
				pdf.savefig()
				plt.close()

	def plotBasis(self, ionName, zeta, n, l, m, axis, minimum=-8, maximum=8, planeValue=0.0, step=0.1, printStatus=False, spectro=True):
		"""Plots cross-section of basis function of ion to pdf.
		All lengths measured in bohr radii (a0).

		Args:
			zeta (int):	Zeta index of Radial
			n (int): Principal quantum number for Radial
			l (int): Orbital angular momentum quantum number for Radial and spherical harmonic
			m (int): Azimuthal quantum number for spherical harmonic
			axis (string) : Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
			minimum (int, optional): Minimum value of coordinates measured in a0; Default is -8
			maximum (int, optional): Maximum value of coordinates measured in a0; Default is +8
			planeValue (double, optional): Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
			step (int, optional): Interval between Cartesian mgrid points, measured in a0; Default is 0.1"""

		ion = self.ions[ionName]
		plotname = 'Basis_'+ionName+'_'+str(zeta)+'_'+str(n)+'_'+str(l)+'_'+str(m)+'_'+axis

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

					# Evaluate value of Radial at mesh point and get psi
					R[i, j] = ion.getRadialValue(zeta, n, l, np.sqrt(space1[i, j]**2 +
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
			ttl = (ionName+' Basis Function for \n \n $\zeta='+str(zeta)+'$, $n='+
				      str(n)+'$, $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+axes[0]+'-'+axes[1]+'$ plane')
			plt.title(ttl)

			# Save to pdf
			pdf.savefig()
			plt.close()
			if printStatus:
				print 'Finished '+plotname+'.pdf'

	@staticmethod
	def plotSPH(cls, l, m, axis, minimum=-8, maximum=8, planeValue=0.0, step=0.1, printStatus=False):
		"""Plots cross-section of spherical harmonic to pdf.
		All lengths measured in bohr radii (a0).

		Input:
		l (int): Orbital angular momentum quantum number for spherical harmonic
		m (int): Azimuthal quantum number for spherical harmonic
		axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int, optional): Minimum value of coordinates measured in a0; Default is -8
		maximum (int, optional): Maximum value of coordinates measured in a0; Default is +8
		planeValue (float, optional): Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
		step (float, optional): Interval between Cartesian mgrid points, measured in a0; Default is 0.1"""

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
			if printStatus:
				print 'Finished '+plotname+'.pdf'

	@staticmethod
	def plotSPH3D(l, m, offset=0.0):

		THETA, PHI = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
		R = np.zeros_like(PHI)

		for i in range(0, len(THETA)):
			for j in range(0, len(PHI)):

				x = np.sin(THETA[i, j]) * np.cos(PHI[i, j])
				y = np.sin(THETA[i, j]) * np.sin(PHI[i, j])
				z = np.cos(THETA[i, j])
				R[i, j] = abs(sph(l, m, x, y, z))


		X = (R + offset) * np.sin(THETA) * np.cos(PHI)
		Y = (R + offset) * np.sin(THETA) * np.sin(PHI)
		Z = (R + offset) * np.cos(THETA)


		fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(12,12))
		im = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, shade=True)


		plt.show()

	def plotPsiCrossSec(self, name, cell, band, axis, minimum, maximum, tolerance=0.0, step=None, planeValue=None, label='', printStatus=False, debug=False):

		plotname = name+'_ChargeDensity_'+axis+'_'+label

		if not step:
			step = cell.gridSpacing

		# Initialise meshes
		# 2D cartesian mesh (x, y, or z axis determined later)
		space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]

		psi2 = np.zeros_like(space1) # Waveunction mesh (psi = R*Y)

		maxPsi2 = 0.0 # Colour plot sets limits to -maxpsi to +maxpsi

		bands = [band]
		if tolerance != 0.0:
			for i in range(0, len(cell.bands)):
				if (abs(cell.bands[band] - cell.bands[i]) < tolerance) and i != band:
					bands.append(i)

		if debug:
			debugString = ''
			for b in bands:
				debugString += str(cell.bands[b])+' '
			print debugString

		# Plot functions to pdf
		with PdfPages('pdfs/' + plotname + '.pdf') as pdf:
			for b in bands:
				# Loop over all mesh points
				for i in range(0, int((maximum - minimum) / step)):
					for j in range(0, int((maximum - minimum) / step)):
						# Use axis variable to determine which axes space1 and space2 refer to
						# Evaluate spherical harmonic at mesh point
						if axis == 'z':
							if not planeValue:
								planeValue = cell.zLength / 2
							psi = cell.givePsi(space2[i, j], space1[i, j], planeValue, band=b)
							plt.xlabel('$x$ / $a_0$')
							plt.ylabel('$y$ / $a_0$')
						if axis == 'y':
							if not planeValue:
								planeValue = cell.yLength / 2
							psi = cell.givePsi(space2[i, j], planeValue, space1[i, j], band=b)
							plt.xlabel('$x$ / $a_0$')
							plt.ylabel('$z$ / $a_0$')
						if axis == 'x':
							if not planeValue:
								planeValue = cell.xLength / 2
							psi = cell.givePsi( planeValue, space2[i, j], space1[i, j], band=b)
							plt.xlabel('$y$ / $a_0$')
							plt.ylabel('$z$ / $a_0$')
					
						psi2[i, j] = psi2[i, j] + abs(psi)**2

						# Update maxpsi
						if abs(psi2[i, j]) > maxPsi2:
							maxPsi2 = psi2[i, j]

			# Setup plot
			plt.imshow(psi2, interpolation='bilinear', origin='center', cmap=plt.cm.Blues,
				          extent=(minimum, maximum, minimum, maximum), vmin=0.0, vmax=maxPsi2)
			plt.colorbar()
			plt.grid()
			axes = ['x', 'y', 'z']
			axes.remove(axis)
			ttl = (cell.name+' Charge Density in $'+axes[0]+'-'+axes[1]+'$ plane')
			plt.title(ttl)

			# Save to pdf
			pdf.savefig()
			plt.close()
			if printStatus:
				print 'Finished '+plotname+'.pdf'

	
