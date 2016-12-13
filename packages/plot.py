
import numpy as np
import datetime as dt

# Import matplotlib packages
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm, colors

from skimage import measure

from sph import sph
from vector import Vector

# Spectroscopic notation dictionary
SPECTRAL = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g',
               5: 'h', 6: 'i', 7: 'j', 8: 'k'}

AXES = ('x', 'y', 'z')


def plot_3d(title, mesh, fraction, x_range, y_range, z_range, step, save_name=None, show=True):
	# Make isosurface at psi2 = fraction * psi2max
	mes = measure.marching_cubes(mesh, fraction*np.max(mesh))
	verts = mes[0]
	faces = mes[1]

	# Set up plot
	fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(10, 10))
	plt.title(title)

	# Plot surface
	ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap=cm.Spectral, lw=0.1)

	# Set axes
	ax.set_xlim3d(x_range[0] / step, x_range[1] / step)
	ax.set_ylim3d(y_range[0] / step, y_range[1] / step)
	ax.set_zlim3d(z_range[0] / step, z_range[1] / step)
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")

	if save_name:
		plt.savefig("figures3D/"+save_name+".png")
	if show:
		plt.show()
	plt.close()


def plot_radials(ions, points=500, printStatus=False, spectro=True):
	"""Plot all radial functions from self.ions to pdf

	Args:
		ions ({string : Ion}): Dict of Ion objects, indexed by ion name
		points (int, opt.): Number of points for plot
		printStatus (bool, opt.): If true, print notification when finishing a plot
		spectro (bool, opt.): If true, use spectroscopic notation
	"""
	timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())

	with PdfPages('pdfs/radials'+timeStamp+'.pdf') as pdf:

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

			for l in ion.radials:
				for zeta in ion.radials[l]:
					# Get Radial and data from ion
					radial = ion.get_radial(l, zeta)
					n = radial.n
					step = radial.cutoff / points
					r = np.arange(0.0, radial.cutoff, step)
					R = np.empty_like(r)
					for i in range(0, len(r)):
						R[i] = radial.getValueCubic(r[i])

					# Add radial info to legend and add to plot
					# If spectro, use spectroscopic notation for legend
					if spectro:
						label = '$\zeta ='+str(zeta)+'$, $'+str(n)+SPECTRAL[l]+'$'
					else:
						label = '$\zeta ='+str(zeta)+'$, $n='+str(n)+'$, $l='+str(l)+'$'
					plt.plot(r, R, label=label)

					if printStatus:
						print "Finished Radial "+ion.ionName+"_"+str(zeta)+"_"+str(n)+"_"+str(l)

			# Add plot to pdf and reset plt
			plt.legend()
			pdf.savefig()
			plt.close()


def plot_sph_2d(l, m, axis, minimum=-8.0, maximum=8.0, planeValue=0.0, step=0.1, printStatus=False):
	"""Plots cross-section of spherical harmonic to pdf.
	All lengths measured in bohr radii (a0).

	Args:
		l (int): Orbital angular momentum quantum number for spherical harmonic
		m (int): Azimuthal quantum number for spherical harmonic
		axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int, opt.): Minimum value of coordinates measured in a0; Default is -8
		maximum (int, opt.): Maximum value of coordinates measured in a0; Default is +8
		planeValue (float, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
		step (float, opt.): Interval between Cartesian mgrid points, measured in a0; Default is 0.1
		printStatus (bool, opt.): If true, print update when plot is finished
	"""

	if axis not in ['x', 'y', 'z']:
		raise ValueError("Axis must be x, y, or z")

	# Initialise meshes
	# 2D cartesian mesh (x, y, or z axis determined later)
	space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]
	Y = np.empty_like(space1, dtype=float)  # Spherical Harmonic mesh

	maxY = 0.1  # Colour plot sets limits to -maxY and +maxY
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
	plotName = 'SPH_'+str(l)+'_'+str(m)+'_'+axis
	with PdfPages('pdfs/'+plotName+'.pdf') as pdf:
		plt.imshow(
			Y, interpolation='bilinear', origin='center', cmap=cm.bwr, extent=(minimum, maximum, minimum, maximum),
			vmin=-maxY, vmax=maxY)
		plt.colorbar()
		plt.grid()
		axes = ['x', 'y', 'z']
		axes.remove(axis)
		ttl = 'Spherical Harmonic for \n \n $l='+str(l)+'$, $m_l='+str(m)+'$ in $'+axes[0]+'-'+axes[1]+'$ plane'
		plt.title(ttl)

		# Save to pdf
		pdf.savefig()
		plt.close()
		if printStatus:
			print 'Finished '+plotName+'.pdf'


def plot_sph_3d(l, m):
	"""Plots 3D spherical harmonic isosurface.

	Args:
		l (int): Degree of spherical harmonic
		m (int): Order of spherical harmonic
	"""

	# Get mesh of angles
	THETA, PHI = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
	# Initialise mesh of 
	SPH = np.zeros_like(PHI, dtype=float)

	# Loop over mesh
	for i in range(0, SPH.shape[0]):
		for j in range(0, SPH.shape[1]):
			# Get cartesian point
			x = np.sin(THETA[i, j]) * np.cos(PHI[i, j])
			y = np.sin(THETA[i, j]) * np.sin(PHI[i, j])
			z = np.cos(THETA[i, j])
			r = Vector(x, y, z)
			SPH[i, j] = abs(sph(l, m, r))
	# Get cartesian mesh
	X = SPH * np.sin(THETA) * np.cos(PHI)
	Y = SPH * np.sin(THETA) * np.sin(PHI)
	Z = SPH * np.cos(THETA)

	# Plot surface
	fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(10, 10))
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")
	ax.set_xlim3d(-1.0, 1.0)
	ax.set_ylim3d(-1.0, 1.0)
	ax.set_zlim3d(-1.0, 1.0)
	title = "Real Spherical Harmonic for Degree $l="+str(l)+"$ and Order $m="+str(m)+"$"
	plt.title(title)
	ax.plot_surface(X, Y, Z, rstride=1, cstride=1)

	plt.show()
	plt.close()


def plot_basis_2d(
		ionName, ion, zeta, n, l, m, axis, minimum=-8, maximum=8, planeValue=0.0, step=0.1, printStatus=False,
		spectro=False):
	"""Plots cross-section of basis function of ion to pdf.
	All lengths measured in bohr radii (a0).

	Args:
		ionName (string): Name of ion file (excluding .ion)
		ion (Ion): Ion object to be plotted
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

	if axis not in ['x', 'y', 'z']:
		raise ValueError("Axis must be x, y, or z")

	# Initialise meshes
	# 2D cartesian mesh (x, y, or z axis determined later)
	space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]

	Y = np.empty_like(space1, dtype=float)  # Spherical Harmonic mesh
	R = np.empty_like(space1, dtype=float)  # Radial Function mesh
	psi = np.empty_like(space1, dtype=float)  # Basis Function mesh (psi = R*Y)

	maxPsi = 0.1  # Colour plot sets limits to -maxPsi to +maxPsi

	# Loop over all mesh points
	for i in range(0, int((maximum - minimum) / step)):
		for j in range(0, int((maximum - minimum) / step)):
			# Use axis variable to determine which axes space1 and space2 refer to
			# Evaluate spherical harmonic at mesh point
			if axis == 'z':
				r = Vector(space2[i, j], space1[i, j], planeValue)
				plt.xlabel('$x$ / $a_0$')
				plt.ylabel('$y$ / $a_0$')
			if axis == 'y':
				r = Vector(space2[i, j], planeValue, space1[i, j])
				plt.xlabel('$x$ / $a_0$')
				plt.ylabel('$z$ / $a_0$')
			if axis == 'x':
				r = Vector(planeValue, space2[i, j], space1[i, j])
				plt.xlabel('$y$ / $a_0$')
				plt.ylabel('$z$ / $a_0$')

			Y[i, j] = sph(l, m, r)
			# Evaluate value of Radial at mesh point and get psi
			distance = abs(r)
			R[i, j] = ion.get_radial_value(l, zeta, distance)
			psi[i, j] = Y[i, j] * R[i, j]

			# Update maxpsi
			if abs(psi[i, j]) > maxPsi:
				maxPsi = abs(psi[i, j])

	if maxPsi == 0:
		raise ValueError("Wavefunction is zero at all points")

	# Setup plot
	timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
	plotName = (
		'Basis_' + ionName + '_' + str(zeta) + '_' + str(n) + '_' + str(l) + '_' + str(m) + '_' + axis + timeStamp)
	with PdfPages('pdfs/' + plotName + '.pdf') as pdf:
		plt.imshow(
			psi, interpolation='bilinear', origin='center', cmap=cm.bwr, extent=(minimum, maximum, minimum, maximum),
			vmin=-maxPsi, vmax=maxPsi)
		plt.colorbar()
		plt.grid()
		axes = ['x', 'y', 'z']
		axes.remove(axis)
		if spectro:
			ttl = (
				ionName+' Basis Function for \n \n $\zeta='+str(zeta)+'$, $'+str(n)+SPECTRAL[l]+'$, $m_l='+str(m)
				+'$ in $'+axes[0]+'-'+axes[1]+'$ plane')
		else:
			ttl = (
				ionName+' Basis Function for \n \n $\zeta='+str(zeta)+'$, $n='+str(n)+'$, $l='+str(l)+'$, $m_l='+str(m)
				+'$ in $'+axes[0]+'-'+axes[1]+'$ plane')
		plt.title(ttl)

		# Save to pdf
		pdf.savefig()
		plt.close()

		if printStatus:
			print 'Finished '+plotName+'.pdf'


def plot_charge_density_gamma_2d(
		cell, E, axis, minimum, maximum, step=None, planeValue=None, label='', printStatus=False):
	"""Plots cross-section of charge density evaluated at gamma-point to pdf.

	All lengths measured in bohr radii (a0).

	Args:
		cell (Cell): Simulation cell to plot
		E (float): Band energy, nearest band energy will be used
		axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int): Minimum value of coordinates
		maximum (int): Maximum value of coordinates
		planeValue (float, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.0
		step (float, opt.): Interval between Cartesian mgrid points, measured in a0;
							Default is cell.gridSpacing
		label (string, opt.): Optional string to append to end of filename
		printStatus (bool, opt.): If true, print update when file is saved
	"""

	if axis not in ['x', 'y', 'z']:
		raise ValueError("Axis must be x, y, or z")

	# If no step given, set to gridSpacing
	if not step:
		step = cell.gridSpacing

	# Initialise meshes
	# 2D cartesian mesh (x, y, or z axis determined later)
	space1, space2 = np.mgrid[minimum:maximum:step, minimum:maximum:step]

	psi2 = np.zeros_like(space1, dtype=float)
	maxPsi2 = 0.0  # Colour plot sets limits to 0 to +maxPsi2

	# Get nearest stored energy at gamma-point to requested energy
	bandEnergy = sorted(cell.get_gamma_energies(), key=lambda t: abs(E - t))[0]

	# Loop over all mesh points
	for i in range(0, int((maximum - minimum) / step)):
		for j in range(0, int((maximum - minimum) / step)):
			# Use axis variable to determine which axes space1 and space2 refer to
			# Evaluate spherical harmonic at mesh point
			if axis == 'z':
				if not planeValue:
					planeValue = cell.zLength / 2
				r = Vector(space2[i, j], space1[i, j], planeValue)
				label1 = '$x$ / $a_0$'
				label2 = '$y$ / $a_0$'
			if axis == 'y':
				if not planeValue:
					planeValue = cell.yLength / 2
				r = Vector(space2[i, j], planeValue, space1[i, j])
				label1 = '$x$ / $a_0$'
				label2 = '$z$ / $a_0$'
			if axis == 'x':
				if not planeValue:
					planeValue = cell.xLength / 2
				r = Vector(planeValue, space2[i, j], space1[i, j])
				label1 = '$y$ / $a_0$'
				label2 = '$z$ / $a_0$'

			psi = cell.get_psi_gamma(bandEnergy, r)
			psi2[i, j] += abs(psi)**2

			# Update maxpsi
			if abs(psi2[i, j]) > maxPsi2:
				maxPsi2 = psi2[i, j]

	if maxPsi2 == 0:
		raise ValueError("Wavefunction is zero at all points")

	# Setup plot
	timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
	plotName = cell.name+'_ChargeDensityGamma_'+axis+'_'+label+timeStamp
	with PdfPages('pdfs/' + plotName + '.pdf') as pdf:
		plt.imshow(
			psi2, interpolation='bilinear', origin='center', cmap=cm.copper,
			extent=(minimum, maximum, minimum, maximum), vmin=0.0, vmax=maxPsi2)
		plt.colorbar()
		plt.xlabel(label1)
		plt.ylabel(label2)
		axes = ['x', 'y', 'z']
		axes.remove(axis)
		ttl = (cell.name+' Charge Density in $'+axes[0]+'-'+axes[1]+'$ plane at $'+axis+'='+str(planeValue)+'$')
		plt.title(ttl)

		# Save to pdf
		pdf.savefig()
		plt.close()
		if printStatus:
			print 'Finished '+plotName+'.pdf'


def plot_charge_density_gamma_3d(
		cell, E, x_range=(0.0, 0.0), y_range=(0.0, 0.0), z_range=(0.0, 0.0), step=0.0, fraction=0.8, alpha=1.0,
		cmap=False, show=True, save=False, debug=False):
	"""Plots charge density isosurface.

	All lengths measured in Bohr radii (a0).

	Args:
		cell (Cell): Simulation cell to plot
		E (float): Band energy
		x_range((float), opt.): Limits of x axis
		y_range((float), opt.): Limits of y axis
		z_range((float), opt.): Limits of z axis
		step (float, opt.): Interval between Cartesian mgrid points; Default is cell.gridSpacing
		fraction (float, opt.): Sets value of isosurface to this fraction of max charge density
		alpha (float, opt.): Transparency of plot surfaces
		cmap (bool, opt.): If true, colour surface opaquely (ignoring alpha) according to z-value
		show (bool, opt.): If true, show plot
		save (bool, opt.): If true, save plot
		debug (bool, opt.): If true, print extra information during runtime
	"""

	# If plot limits not given, set to limits of cell
	if x_range == (0.0, 0.0):
		x_range = (0.0, cell.vector.x)
	if y_range == (0.0, 0.0):
		y_range = (0.0, cell.vector.y)
	if z_range == (0.0, 0.0):
		z_range = (0.0, cell.vector.z)

	# If step not given, set to cell.gridSpacing
	if step == 0.0:
		step = cell.grid_spacing

	# Get nearest stored energy at gamma-point to requested energy
	bandEnergy = sorted(cell.get_gamma_energies(), key=lambda t: abs(E - t))[0]

	psi = cell.get_psi_grid(Vector.zero(), bandEnergy, debug=debug)
	psi2 = abs(psi)**2
	max_psi2 = np.max(psi2)
	if max_psi2 == 0.0:
		raise ValueError("Wavefunction is zero at all points")

	title = (
		cell.name+' Charge Density Isosurface at '+str(fraction)+' of Maximum Density at \n Gamma Point and Energy '
		+str(bandEnergy)+' eV relative to the Fermi Level')

	# Save plot as png
	if save:
		timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
		save_name = cell.name+"_ChargeDensity3D_"+str(fraction)+"_"+str(bandEnergy)+"_"+timeStamp
	else:
		save_name = None

	plot_3d(title, psi2, fraction, x_range, y_range, z_range, step, save_name=save_name, show=show)


def plot_ldos_2d(
		cell, min_E, max_E, T, axis, planeValue, minimum, maximum, step=None, interpolation='cubic',
		printStatus=False, debug=False):
	"""Plots cross-section of charge density to pdf.

	All lengths measured in bohr radii (a0).

	Args:
		cell (Cell): Simulation cell to plot
		min_E (float): Minimum energy
		max_E (float): Maximum energy
		T (float): Absolute temperature in K
		axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int): Minimum value of coordinates
		maximum (int): Maximum value of coordinates
		planeValue (float, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.0
		step (float, opt.): Interval between Cartesian mgrid points, measured in a0;
							Default is cell.gridSpacing
		interpolation (string, opt.): Method of interpolation; possible arguments are 'cubic' (default) and 'linear'
		printStatus (bool, opt.): If true, print update when file is saved
		debug (bool, opt.): If true, print extra information during runtime
	"""

	if axis not in AXES:
		raise ValueError("Axis must be x, y, or z")

	# If no step given, set to gridSpacing
	if not step:
		step = cell.grid_spacing

	ldos_3d = cell.get_ldos_grid(min_E, max_E, T, debug=debug)

	if np.max(ldos_3d) == 0.0:
		raise ValueError("LDoS is zero at all points")


	timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
	save_name = cell.name + '_LDoS2D_' + axis + '_' + timeStamp

	axes = ['x', 'y', 'z']
	axes.remove(axis)
	title = cell.name+' LDoS in $'+axes[0]+'-'+axes[1]+'$ plane at $'+axis+'='+str(planeValue)+'a_0$'

	plot_2d(cell, ldos_3d, title, save_name, axis, planeValue, minimum, maximum)

	if printStatus:
		print 'Finished ' + save_name + '.pdf'


def plot_2d(cell, mesh_3d, title, save_name, axis, plane_value, minimum, maximum):

	if axis == 'x':
		index = np.where(cell.xMesh == plane_value)[0][0]
		mesh_2d = mesh_3d[index, :, :]
		label1 = 'z'
		label2 = 'y'
	elif axis == 'y':
		index = np.where(cell.yMesh == plane_value)[1][0]
		mesh_2d = mesh_3d[:, index, :]
		label1 = 'z'
		label2 = 'x'
	elif axis == 'z':
		index = np.where(cell.zMesh == plane_value)[2][0]
		mesh_2d = mesh_3d[:, :, index]
		label1 = 'y'
		label2 = 'x'

	with PdfPages('figures2D/'+save_name+'.pdf') as pdf:
		plt.imshow(
				mesh_2d, interpolation='bilinear', origin='center', cmap=cm.copper,
				extent=(minimum, maximum, minimum, maximum))
		plt.colorbar()
		plt.title(title)
		plt.xlabel(label1)
		plt.ylabel(label2)

		# Save to pdf
		pdf.savefig()
		plt.close()


def plot_ldos_3d(
		cell, min_E, max_E, T, x_range=(0.0, 0.0), y_range=(0.0, 0.0), z_range=(0.0, 0.0), step=0.0, fraction=0.8, alpha=1.0,
		show=True, save=False, debug=False, recalculate=False):
	"""Plots charge density isosurface.

	All lengths measured in Bohr radii (a0).

	Args:
		cell (Cell): Simulation cell to plot
		min_E (float): Minimum of energy range
		max_E (float): Maximum of energy range
		x_range((float), opt.): Limits of x axis
		y_range((float), opt.): Limits of y axis
		z_range((float), opt.): Limits of z axis
		step (float, opt.): Interval between Cartesian mgrid points; Default is cell.gridSpacing
		fraction (float, opt.): Sets value of isosurface to this fraction of max charge density
		alpha (float, opt.): Transparency of plot surfaces
		show (bool, opt.): If true, show plot
		save (bool, opt.): If true, save plot
		debug (bool, opt.): If true, print extra information during runtime
	"""

	# If plot limits not given, set to limits of cell
	if x_range == (0.0, 0.0):
		x_range = (0.0, cell.vector.x)
	if y_range == (0.0, 0.0):
		y_range = (0.0, cell.vector.y)
	if z_range == (0.0, 0.0):
		z_range = (0.0, cell.vector.z)

	# If step not given, set to cell.grid_spacing
	if step == 0.0:
		step = cell.grid_spacing

	# Get energy relative to Fermi level
	min_EAbsolute = min_E + cell.fermiLevel
	max_EAbsolute = max_E + cell.fermiLevel

	# Cartesian mesh
	ldos = cell.get_ldos_grid(min_EAbsolute, max_EAbsolute, T, debug=debug, recalculate=recalculate)
	max_ldos = np.max(ldos)

	if max_ldos == 0.0:
		raise ValueError("LDoS is zero at all points")

	title = (
		cell.name+' LDoS Isosurface at '+str(fraction)+' of Maximum Density for \n Energies from '
		+str(min_E)+' eV to '+str(max_E)+' eV relative to the Fermi Level')

	# Save plot as png
	if save:
		timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
		save_name = cell.name+"_LDoS3D_"+str(fraction)+"_"+str(min_EAbsolute)+"_"+str(max_EAbsolute)+"_"+timeStamp
	else:
		save_name = None

	plot_3d(title, ldos, fraction, x_range, y_range, z_range, step, save_name=save_name, show=show)
