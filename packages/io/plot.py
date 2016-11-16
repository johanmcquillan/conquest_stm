import numpy as np

# Import matplotlib packages
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm, colors

from skimage import measure

from packages.sph import sph
from packages.smartDict import SmartDict

# Spectroscopic notation dictionary
SPECTRAL = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g',
               5: 'h', 6: 'i', 7: 'j', 8: 'k'}

def plotRadials(ions, points=500, printStatus=False, spectro=True):
	"""Plot all radial functions from self.ions to 'self.filename'_radials.pdf

	Args:
		points (int, opt.): Number of points for plot
		printStatus (boolean, opt.): If true, print notification when finishing a plot
		spectro (boolean, opt.): If true, use spectroscopic notation
	"""
	with PdfPages('pdfs/radials.pdf') as pdf:

		# Plot all functions for the same ion on one graph
		ionNames = sorted(ions.keys())
		for ionName in ionNames:
			ion = ions[ionName]

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
							label = '$\zeta ='+str(zeta)+'$, $'+str(n)+self.SPECTRAL[l]+'$'
						else:
							label = '$\zeta ='+str(zeta)+'$, $n='+str(n)+'$, $l='+str(l)+'$'
						plt.plot(r, R, label=label)

						if printStatus:
							print "Finished Radial "+ion.ionName+"_"+str(zeta)+"_"+str(n)+"_"+str(l)

			# Add plot to pdf and reset plt
			plt.legend()
			pdf.savefig()
			plt.close()

def plotSPH2D(l, m, axis, minimum=-8.0, maximum=8.0, planeValue=0.0, step=0.1, printStatus=False):
	"""Plots cross-section of spherical harmonic to pdf.
	All lengths measured in bohr radii (a0).

	Args:
		l (int): Orbital angular momentum quantum number for spherical harmonic
		m (int): Azimuthal quantum number for spherical harmonic
		axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int, optional): Minimum value of coordinates measured in a0; Default is -8
		maximum (int, optional): Maximum value of coordinates measured in a0; Default is +8
		planeValue (float, optional): Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
		step (float, optional): Interval between Cartesian mgrid points, measured in a0; Default is 0.1
	"""

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
		plt.imshow(Y, interpolation='bilinear', origin='center', cmap=cm.bwr,
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

def plotSPH3D(l, m):
	"""Plots 3D spherical harmonic isosurface.

	Input:
		l (int): Degree of spherical harmonic
		m (int): Order of spherical harmonic
	"""

	# Get mesh of angles
	THETA, PHI = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
	# Initialise mesh of 
	SPH = np.zeros_like(PHI)

	# Loop over mesh
	for i in range(0, SPH.shape[0]):
		for j in range(0, SPH.shape[1]):
			# Get cartesian point
			x = np.sin(THETA[i, j]) * np.cos(PHI[i, j])
			y = np.sin(THETA[i, j]) * np.sin(PHI[i, j])
			z = np.cos(THETA[i, j])
			SPH[i, j] = abs(sph(l, m, x, y, z))

	# Get cartesian mesh
	X = SPH * np.sin(THETA) * np.cos(PHI)
	Y = SPH * np.sin(THETA) * np.sin(PHI)
	Z = SPH * np.cos(THETA)

	# Plot surface
	fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(10,10))
	ax.set_xlim3d(-1.0, 1.0)
	ax.set_ylim3d(-1.0, 1.0)
	ax.set_zlim3d(-1.0, 1.0)
	im = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, shade=False)

	plt.show()
	plt.close()

def plotBasis2D(ionName, ion, zeta, n, l, m, axis, minimum=-8, maximum=8, planeValue=0.0, step=0.1, printStatus=False, spectro=True):
	"""Plots cross-section of basis function of ion to pdf.
	All lengths measured in bohr radii (a0).

	Args:
		zeta (int):	Zeta index of Radial
		n (int): Principal quantum number for Radial
		l (int): Orbital angular momentum quantum number for Radial and spherical harmonic
		m (int): Azimuthal quantum number for spherical harmonic
		axis (string) : Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int, opt.): Minimum value of coordinates measured in a0; Default is -8
		maximum (int, opt.): Maximum value of coordinates measured in a0; Default is +8
		planeValue (double, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
		step (int, opt.): Interval between Cartesian mgrid points, measured in a0; Default is 0.1
		printStatus (boolean, opt.): If true, print notification when finishing a plot
		spectro (boolean, opt.): If true, use spectroscopic notation
	"""

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
				distance = np.sqrt(space1[i, j]**2 + space2[i, j]**2 + planeValue**2)
				R[i, j] = ion.getRadialValue(zeta, n, l, distance)
				psi[i, j] = Y[i, j] * R[i, j]

				# Update maxpsi
				if abs(psi[i, j]) > maxpsi:
					maxpsi = abs(psi[i, j])

		# Setup plot
		plt.imshow(psi, interpolation='bilinear', origin='center', cmap=cm.bwr,
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

def plotChargeDensity2D(cell, bandNumber, axis, minimum, maximum, step=None,
	                       planeValue=None, normalise=False, label='', printStatus=False, debug=False):
	"""Plots cross-section of charge density to pdf.
	All lengths measured in bohr radii (a0).

	Args:
		cell (Cell): Simulation cell to plot
		bandNumber (int): Band number to plot
		axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int): Minimum value of coordinates
		maximum (int): Maximum value of coordinates
		planeValue (float, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.0
		normalise (bool, opt.): If true, normalise basis coefficients before plot
		step (float, opt.): Interval between Cartesian mgrid points, measured in a0;
							Default is cell.gridSpacing
		label (string, opt.): Optional string to append to end of filename
		debug (bool, opt.): If true, print extra information during runtime
	"""

	# If no step given, set to gridSpacing
	if not step:
		step = cell.gridSpacing
	# Normalise basis coefficients
	if normalise:
		cell.normaliseBand(bandNumber)

	plotname = cell.name+'_ChargeDensity_'+axis+'_'+label

	# Initialise meshes
	# 2D cartesian mesh (x, y, or z axis determined later)
	space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]

	psi2 = np.zeros_like(space1) # Waveunction mesh (psi = R*Y)

	maxPsi2 = 0.0 # Colour plot sets limits to -maxpsi to +maxpsi

	# Print debug info
	if debug:
		debugString = ''
		for b in bandNumbers:
			debugString += str(cell.bands[b])+' '
		print debugString

	# Plot functions to pdf
	with PdfPages('pdfs/' + plotname + '.pdf') as pdf:
		# Loop over all mesh points
		for i in range(0, int((maximum - minimum) / step)):
			for j in range(0, int((maximum - minimum) / step)):
				# Use axis variable to determine which axes space1 and space2 refer to
				# Evaluate spherical harmonic at mesh point
				if axis == 'z':
					if not planeValue:
						planeValue = cell.zLength / 2
					psi = cell.givePsi(space2[i, j], space1[i, j], planeValue, bandNumber=bandNumber)
					plt.xlabel('$x$ / $a_0$')
					plt.ylabel('$y$ / $a_0$')
				if axis == 'y':
					if not planeValue:
						planeValue = cell.yLength / 2
					psi = cell.givePsi(space2[i, j], planeValue, space1[i, j], bandNumber=bandNumber)
					plt.xlabel('$x$ / $a_0$')
					plt.ylabel('$z$ / $a_0$')
				if axis == 'x':
					if not planeValue:
						planeValue = cell.xLength / 2
					psi = cell.givePsi( planeValue, space2[i, j], space1[i, j], bandNumber=bandNumber)
					plt.xlabel('$y$ / $a_0$')
					plt.ylabel('$z$ / $a_0$')
			
				psi2[i, j] = psi2[i, j] + abs(psi)**2

				# Update maxpsi
				if abs(psi2[i, j]) > maxPsi2:
					maxPsi2 = psi2[i, j]

		# Setup plot
		plt.imshow(psi2, interpolation='bilinear', origin='center', cmap=cm.Blues,
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

def plotChargeDensity3D(cell, bandNumber, xrange=(0.0, 0.0), yrange=(0.0, 0.0), zrange=(0.0, 0.0),
	                       step=0.0, fraction=0.8, alpha=1.0, cmap=False, normalise=False, show=True, save=False):
	"""Plots charge density isosurface.

	All lengths measured in Bohr radii (a0).

	Args:
		cell (Cell): Simulation cell to plot
		bandNumber (int): Band number to plot
		xrange((float), opt.): Limits of x axis
		yrange((float), opt.): Limits of y axis
		zrange((float), opt.): Limits of z axis
		step (float, opt.): Interval between Cartesian mgrid points; Default is cell.gridSpacing
		fraction (float, opt.): Sets value of isosurface to this fraction of max charge density
		alpha (float, opt.): Transparency of plot surfaces
		cmap (bool, opt.): If true, colour surface opaquely (ignoring alpha) according to z-value
		normalise (bool, opt.): If true, normalise basis coefficients before plot
		show (bool, opt.): If true, show plot
		save (bool, opt.): If true, save plot
	"""

	bandEnergy = cell.bands[bandNumber]
	# If plot limits not given, set to limits of cell
	if xrange == (0.0, 0.0):
		xrange = (0.0, cell.xLength)
	if yrange == (0.0, 0.0):
		yrange = (0.0, cell.yLength)
	if zrange == (0.0, 0.0):
		zrange = (0.0, cell.zLength)

	# If step not given, set to cell.gridSpacing
	if step == 0.0:
		step = cell.gridSpacing

	if normalise:
		cell.normaliseBand(bandNumber)

	# Cartesian mesh
	X, Y, Z = np.mgrid[xrange[0]:xrange[1]:step,
	                      yrange[0]:yrange[1]:step,
	                      zrange[0]:zrange[1]:step]

	psi2 = np.zeros_like(X, dtype=float)
	psi2max = 0.0

	# Loop over all mesh points
	for i in range(int((xrange[1] - xrange[0]) / step)):
		for j in range(int((yrange[1] - yrange[0]) / step)):
			for k in range(int((zrange[1] - zrange[0]) / step)):
				# Get coordinates
				x = X[i, j, k]
				y = Y[i, j, k]
				z = Z[i, j, k]

				# Calculate wavefunction
				psi = cell.givePsi(x, y, z, bandNumber=bandNumber)

				# Get charge density
				psi2[i, j, k] = abs(psi)**2

				# Set max value
				if psi2[i, j, k] > psi2max:
					psi2max = psi2[i, j, k]

	# Make isosurface at psi2 = fraction * psi2max
	mes = measure.marching_cubes(psi2, fraction*psi2max)
	verts = mes[0]
	faces = mes[1]

	# Set up plot
	fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(10,10))
	title = (cell.name+' Charge Density Isosurface at '+str(fraction)+' of Maximum Density for \n Band Energy '+str(bandEnergy)+' Ha with Fermi Level '+str(cell.fermiLevel)+' Ha')
	plt.title(title)

	# Set axes
	ax.set_xlim3d(xrange[0]/step, xrange[1]/step)
	ax.set_ylim3d(yrange[0]/step, yrange[1]/step)
	ax.set_zlim3d(zrange[0]/step, zrange[1]/step)
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")

	# Plot surface
	if cmap:
		ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap=cm.Spectral, lw=0.1)
	else:
		ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], color=(1,0,0,alpha), lw=0.1)

	# Save plot as png
	if save:
		saveName = cell.name+"_ChargeDensity3D_"+str(fraction)+"_"+str(bandEnergy)
		plt.savefig("figures3D/"+saveName+".png")
	# Show plot
	if show:
		plt.show()
	plt.close()
