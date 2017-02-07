
import numpy as np
import datetime as dt

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib import colors as clrs

from skimage import measure

from sph import sph
from vector import Vector

# Spectroscopic notation dictionary
SPECTRAL = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g',
               5: 'h', 6: 'i', 7: 'j', 8: 'k'}

AXES = ('x', 'y', 'z')


def plot_2d(cell, mesh_3d, title, save_name, axis, plane_value):

	if axis == 'x':
		index = np.where(cell.x_mesh == plane_value)[0][0]
		mesh_2d = mesh_3d[index, :, :]
		label1 = 'z'
		label2 = 'y'
	elif axis == 'y':
		index = np.where(cell.y_mesh == plane_value)[1][0]
		mesh_2d = mesh_3d[:, index, :]
		label1 = 'z'
		label2 = 'x'
	elif axis == 'z':
		index = np.where(cell.z_mesh == plane_value)[2][0]
		mesh_2d = mesh_3d[:, :, index]
		label1 = 'y'
		label2 = 'x'
	else:
		raise ValueError("Axis must be x, y, or z")

	with PdfPages('figures2D/'+save_name+'.pdf') as pdf:
		plt.imshow(
				mesh_2d, interpolation='bilinear', origin='lower', cmap=cm.copper,
				extent=(0, mesh_2d.shape[0], 0, mesh_2d.shape[1]))
		plt.colorbar()
		plt.title(title)
		plt.xlabel(label1)
		plt.ylabel(label2)

		# Save to pdf
		pdf.savefig()
		plt.close()


def plot_3d_scatter_mask(title, mesh, delta=0.0, zero_surf=False):

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	plt.title(title)

	points = []

	for ijk in np.ndindex(mesh.shape):
		# print ijk, mesh[ijk], (type(mesh) is not np.ma.MaskedArray or mesh[ijk] is not np.ma.masked), (delta == 0.0 or abs(mesh[ijk]) <= delta), (zero_surf or mesh[ijk] != 0)
		if (type(mesh) is not np.ma.MaskedArray or mesh[ijk] is not np.ma.masked) and (delta == 0.0 or abs(mesh[ijk]) <= delta) and (zero_surf or mesh[ijk] != 0):
			i, j, k = ijk
			points.append([i, j, k, mesh[ijk]])
			print (i, j, k), mesh[ijk]
			#print i, j, k, mesh[ijk], delta

	np_points = np.transpose(np.array(points))

	if len(points) == 0 or np.max(np_points) == 0:
		raise ValueError("No Points!")
	else:
		print np_points.shape, np.max(np_points[3])

	#colmap = cm.ScalarMappable()
	#colmap.set_array(np_points[3])

	ax.scatter(np_points[0], np_points[1], np_points[2], c=cm.Reds(abs(np_points[3])/np.max(np_points[3])))

	plt.show()


def plot_3d_isosurface(title, mesh, fraction, x_range, y_range, z_range, step, save_name=None, show=True, top_down=False):
	"""Plot isosurface of 3D mesh.

	Args:
		fraction (float, opt.): Sets value of isosurface to this fraction of max charge density
		show (bool, opt.): If true, show plot
		save (bool, opt.): If true, save plot
		debug (bool, opt.): If true, print extra information during runtime
	"""

	mes = measure.marching_cubes(mesh, fraction*np.max(mesh))
	verts = mes[0]
	faces = mes[1]

	# Set up plot
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	plt.title(title)

	ax.set_xlim3d(x_range[0] / step, x_range[1] / step)
	ax.set_ylim3d(y_range[0] / step, y_range[1] / step)
	ax.set_zlim3d(z_range[0] / step, z_range[1] / step)

	if top_down:
		ax.view_init(elev=90, azim=-90)
		ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap=cm.Greys_r, antialiased=False, lw=0.0, vmin=55)
	else:
		# Set axes
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_zlabel("z")
		ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap=cm.Spectral, antialiased=True, lw=0.1)

	if save_name:
		plt.savefig("figures3D/"+save_name+".png")
	if show:
		plt.show()
	plt.close()


def plot_vector_field(cell, vector_field):

	# fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

	fig = plt.figure()
	ax = fig.gca(projection='3d')

	# for i in range(vector_field.shape[0]):
	# 	for j in range(vector_field.shape[1]):
	# 		for k in range(vector_field.shape[2]):
	# 			print vector_field[i, j, k]

	ax.quiver(cell.x_mesh, cell.y_mesh, cell.z_mesh, vector_field[0], vector_field[1], vector_field[2])

	plt.show()


def plot_radials(ions, points=500, print_status=False, spectro=True, interpolation='cubic'):
	"""Plot all radial functions from self.ions to pdf

	Args:
		ions ({string : Ion}): Dict of Ion objects, indexed by ion name
		points (int, opt.): Number of points for plot
		print_status (bool, opt.): If true, print notification when finishing a plot
		spectro (bool, opt.): If true, use spectroscopic notation
	"""
	timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())

	with PdfPages('pdfs/radials'+timeStamp+'.pdf') as pdf:

		# Plot all functions for the same ion on one graph
		ion_names = sorted(ions.keys())
		for ion_name in ion_names:
			ion = ions[ion_name]

			# Setup plot
			fig, ax = plt.subplots(figsize=(10, 5))

			plt.title('Radial Functions for '+ion_name+'.ion')
			ax.set_xlabel('Radial Distance, $r$ / $a_0$')
			ax.set_ylabel('$R_{\zeta l}(r)$')
			plt.minorticks_on()
			plt.grid(b=True, which='major', alpha=0.25, linestyle='-')
			plt.grid(b=True, which='minor', alpha=0.05, linestyle='-')

			# Loop over all radial functions
			maxy = 0
			my = 0
			for l in ion.radials:
				for zeta in ion.radials[l]:
					# Get Radial and data from ion
					radial = ion.get_radial(l, zeta)
					my = max(radial.radial_function_values)
					if my > maxy:
						maxy = my
					n = radial.n
					step = radial.cutoff / points
					r = np.arange(0.0, radial.cutoff, step)
					R = np.empty_like(r)
					for i in range(0, len(r)):
						R[i] = radial.get_value(r[i], interpolation=interpolation)

					# Add radial info to legend and add to plot
					# If spectro, use spectroscopic notation for legend
					if spectro:
						label = '$\zeta ='+str(zeta)+'$; $'+str(n)+SPECTRAL[l]+'$'
					else:
						label = '$\zeta ='+str(zeta)+'$; $n='+str(n)+'$, $l='+str(l)+'$'
					ax.plot(r, R, label=label)

					if print_status:
						print "Finished Radial "+ion.ion_name+"_"+str(zeta)+"_"+str(n)+"_"+str(l)
			ymax = 0.2 * int(maxy / 0.2 + 2)
			ax.set_ylim(ymin=0, ymax=ymax)

			# Add plot to pdf and reset plt
			plt.legend()
			pdf.savefig()
			plt.close()


def plot_sph_2d(l, m, axis, minimum=-8.0, maximum=8.0, plane_value=0.0, step=0.1, print_status=False):
	"""Plots cross-section of spherical harmonic to pdf.
	All lengths measured in bohr radii (a0).

	Args:
		l (int): Orbital angular momentum quantum number for spherical harmonic
		m (int): Azimuthal quantum number for spherical harmonic
		axis (string): Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int, opt.): Minimum value of coordinates measured in a0; Default is -8
		maximum (int, opt.): Maximum value of coordinates measured in a0; Default is +8
		plane_value (float, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
		step (float, opt.): Interval between Cartesian mgrid points, measured in a0; Default is 0.1
		print_status (bool, opt.): If true, print update when plot is finished
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
				r = Vector(space2[i, j], space1[i, j], plane_value)
				plt.xlabel('$x$ / $a_0$')
				plt.ylabel('$y$ / $a_0$')
			if axis == 'y':
				r = Vector(space2[i, j], plane_value, space1[i, j])
				plt.xlabel('$x$ / $a_0$')
				plt.ylabel('$z$ / $a_0$')
			if axis == 'x':
				r = Vector(plane_value, space2[i, j], space1[i, j])
				plt.xlabel('$y$ / $a_0$')
				plt.ylabel('$z$ / $a_0$')

			Y[i, j] = sph(l, m, r)
			# Update maxY
			if abs(Y[i, j]) > maxY:
				maxY = abs(Y[i, j])

	# Setup plot
	plot_name = 'SPH_'+str(l)+'_'+str(m)+'_'+axis
	with PdfPages('pdfs/'+plot_name+'.pdf') as pdf:
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
		if print_status:
			print 'Finished '+plot_name+'.pdf'


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
		ion_name, ion, zeta, n, l, m, axis, minimum=-8, maximum=8, plane_value=0.0, step=0.1, print_status=False,
		spectro=False):
	"""Plots cross-section of basis function of ion to pdf.
	All lengths measured in bohr radii (a0).

	Args:
		ion_name (string): Name of ion file (excluding .ion)
		ion (Ion): Ion object to be plotted
		zeta (int):	Zeta index of Radial
		n (int): Principal quantum number for Radial
		l (int): Orbital angular momentum quantum number for Radial and spherical harmonic
		m (int): Azimuthal quantum number for spherical harmonic
		axis (string) : Cartesian axis ('x', 'y', or 'z') to set to constant value given by planeValue
		minimum (int, opt.): Minimum value of coordinates measured in a0; Default is -8
		maximum (int, opt.): Maximum value of coordinates measured in a0; Default is +8
		plane_value (double, opt.): Constant value assigned to Cartesian coordinate given by axis; Default is 0.00001
		step (int, opt.): Interval between Cartesian mgrid points, measured in a0; Default is 0.1
		print_status (boolean, opt.): If true, print notification when finishing a plot
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

	max_psi = 0.1  # Colour plot sets limits to -max_psi to +max_psi

	# Loop over all mesh points
	for i in range(0, int((maximum - minimum) / step)):
		for j in range(0, int((maximum - minimum) / step)):
			# Use axis variable to determine which axes space1 and space2 refer to
			# Evaluate spherical harmonic at mesh point
			if axis == 'z':
				r = Vector(space2[i, j], space1[i, j], plane_value)
				plt.xlabel('$x$ / $a_0$')
				plt.ylabel('$y$ / $a_0$')
			if axis == 'y':
				r = Vector(space2[i, j], plane_value, space1[i, j])
				plt.xlabel('$x$ / $a_0$')
				plt.ylabel('$z$ / $a_0$')
			if axis == 'x':
				r = Vector(plane_value, space2[i, j], space1[i, j])
				plt.xlabel('$y$ / $a_0$')
				plt.ylabel('$z$ / $a_0$')

			Y[i, j] = sph(l, m, r)
			# Evaluate value of Radial at mesh point and get psi
			distance = abs(r)
			R[i, j] = ion.get_radial_value(l, zeta, distance)
			psi[i, j] = Y[i, j] * R[i, j]

			# Update maxpsi
			if abs(psi[i, j]) > max_psi:
				max_psi = abs(psi[i, j])

	if max_psi == 0:
		raise ValueError("Wavefunction is zero at all points")

	# Setup plot
	time_stamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
	plot_name = (
		'Basis_' + ion_name + '_' + str(zeta) + '_' + str(n) + '_' + str(l) + '_' + str(m) + '_' + axis + time_stamp)
	with PdfPages('pdfs/' + plot_name + '.pdf') as pdf:
		plt.imshow(
			psi, interpolation='bilinear', origin='center', cmap=cm.bwr, extent=(minimum, maximum, minimum, maximum),
			vmin=-max_psi, vmax=max_psi)
		plt.colorbar()
		plt.grid()
		axes = ['x', 'y', 'z']
		axes.remove(axis)
		if spectro:
			ttl = (
				ion_name + ' Basis Function for \n \n $\zeta=' + str(zeta) + '$, $' + str(n) + SPECTRAL[l] + '$, $m_l=' + str(m)
				+ '$ in $' + axes[0] + '-' + axes[1] + '$ plane')
		else:
			ttl = (
				ion_name + ' Basis Function for \n \n $\zeta=' + str(zeta) + '$, $n=' + str(n) + '$, $l=' + str(l) + '$, $m_l=' + str(m)
				+ '$ in $' + axes[0] + '-' + axes[1] + '$ plane')
		plt.title(ttl)

		# Save to pdf
		pdf.savefig()
		plt.close()

		if print_status:
			print 'Finished '+plot_name+'.pdf'


def plot_ldos_2d(
		cell, min_E, max_E, T, axis, planeValue, minimum, maximum, step=None, interpolation='cubic',
		print_status=False, debug=False):
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
		print_status (bool, opt.): If true, print update when file is saved
		debug (bool, opt.): If true, print extra information during runtime
	"""

	if axis not in AXES:
		raise ValueError("Axis must be x, y, or z")

	# If no step given, set to gridSpacing
	if not step:
		step = cell.grid_spacing

	# Get absolute energy
	min_E_abs = min_E + cell.fermiLevel
	max_E_abs = max_E + cell.fermiLevel

	ldos_3d = cell.get_ldos_grid(min_E_abs, max_E_abs, T, debug=debug)

	if np.max(ldos_3d) == 0.0:
		raise ValueError("LDoS is zero at all points")
	else:
		index = np.where(ldos_3d == 16.9705557301)
		i = index[0]
		j = index[1]
		k = index[2]

	timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
	save_name = cell.name + '_LDoS2D_' + axis + '_' + timeStamp

	axes = ['x', 'y', 'z']
	axes.remove(axis)
	title = cell.name+' LDoS in $'+axes[0]+'-'+axes[1]+'$ plane at $'+axis+'='+str(planeValue)+'a_0$'

	plot_2d(cell, ldos_3d, title, save_name, axis, planeValue, minimum, maximum)

	if print_status:
		print 'Finished ' + save_name + '.pdf'


def plot_ldos_3d(
		cell, min_E, max_E, T, step=0.0, fraction=0.8,
		show=True, save=False, recalculate=False, vectorised=True, top_down=False, debug=False):
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
	x_range = (0.0, cell.vector.x)
	y_range = (0.0, cell.vector.y)
	z_range = (0.0, cell.vector.z)

	# If step not given, set to cell.grid_spacing
	if step == 0.0:
		step = cell.grid_spacing

	# Get absolute energy
	min_E_abs = min_E + cell.fermi_level
	max_E_abs = max_E + cell.fermi_level

	# Cartesian mesh
	ldos = cell.get_ldos_grid(min_E_abs, max_E_abs, T, vectorised=vectorised, debug=debug, recalculate=recalculate)
	max_ldos = np.max(ldos)

	if max_ldos == 0.0:
		raise ValueError("LDoS is zero at all points")

	title = (
		cell.name+' LDoS Isosurface at '+str(fraction)+' of Maximum Density for \n Energies from '
		+str(min_E)+' eV to '+str(max_E)+' eV relative to the Fermi Level')

	# Save plot as png
	if save:
		timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
		save_name = cell.name+"_LDoS3D_"+str(fraction)+"_"+str(min_E_abs)+"_"+str(max_E_abs)+"_"+timeStamp
	else:
		save_name = None

	plot_3d_isosurface(title, ldos, fraction, x_range, y_range, z_range, step, save_name=save_name, show=show, top_down=top_down)

def plot_current_2d_iso(
		cell, z, V, T, tip_work_func, tip_energy, delta_s, interpolation='cubic',
		print_status=False, recalculate=False, show=True, partial_surface=False, debug=False):
	"""Plots cross-section of charge density to pdf.

	All lengths measured in bohr radii (a0).

	Args:
		cell (Cell): Simulation cell to plot
		min_E (float): Minimum energy
		max_E (float): Maximum energy
		T (float): Absolute temperature in K
		interpolation (string, opt.): Method of interpolation; possible arguments are 'cubic' (default) and 'linear'
		print_status (bool, opt.): If true, print update when file is saved
		debug (bool, opt.): If true, print extra information during runtime
	"""
	current = cell.get_current_scan_iso(z, V, T, tip_work_func, tip_energy, delta_s, recalculate=recalculate, debug=debug, partial_surface=partial_surface)

	timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
	save_name = cell.name + '_current_' + str(z) +'_' + str(V) + '_' + str(T) + timeStamp

	title = cell.name+' STM scan at $V={:.2}V$ at $z={}a_0$'.format(V, z)

	with PdfPages('figures2D/'+save_name+'.pdf') as pdf:
		plt.imshow(current, interpolation='bilinear', origin='lower', cmap=cm.afmhot)
		plt.colorbar()
		plt.title(title)
		plt.xlabel('y')
		plt.ylabel('x')

		# Save to pdf
		pdf.savefig()

		if show:
			plt.show()

		plt.close()

	if print_status:
		print 'Finished ' + save_name + '.pdf'
		
		
def plot_current_2d_plane(
		cell, z, wf_height, V, T, tip_work_func, tip_energy, delta_s, interpolation='cubic',
		print_status=False, recalculate=False, show=True, debug=False):
	"""Plots cross-section of charge density to pdf.

	All lengths measured in bohr radii (a0).

	Args:
		cell (Cell): Simulation cell to plot
		T (float): Absolute temperature in K
		interpolation (string, opt.): Method of interpolation; possible arguments are 'cubic' (default) and 'linear'
		print_status (bool, opt.): If true, print update when file is saved
		debug (bool, opt.): If true, print extra information during runtime
	"""
	current = cell.get_current_scan_plane(z, wf_height, V, T, tip_work_func, tip_energy, recalculate=recalculate, debug=debug)
	
	timeStamp = '_{:%Y-%m-%d-%H-%M-%S}'.format(dt.datetime.now())
	save_name = cell.name + '_current_' + str(z) +'_' + str(V) + '_' + str(T) + timeStamp

	title = cell.name+r' STM scan at $V={:.2}V$ at $z={}a_0$, with $\psi$ integrated at $z={}$'.format(V, z, wf_height)

	with PdfPages('figures2D/'+save_name+'.pdf') as pdf:
		plt.imshow(current, interpolation='bilinear', origin='lower', cmap=cm.copper)
		plt.colorbar()
		plt.title(title)
		plt.xlabel('y')
		plt.ylabel('x')

		# Save to pdf
		pdf.savefig()

		if show:
			plt.show()

		plt.close()

	if print_status:
		print 'Finished ' + save_name + '.pdf'
