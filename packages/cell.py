
import os
import numpy as np
import math
from skimage import measure

from sph import sph
from smart_dict import SmartDict
from vector import Vector, KVector
from safe_open import safe_open

BOLTZMANN = 8.6173303E-5  # Boltzmann's Constant in eV/K
ELECTRON_MASS = 9.10938E-31  # Electron Mass in kg
H_BAR = 4.135667662E-15  # Reduced Planck's Constant in eV.s

MESH_FOLDER = "temp/"
SUPPORT_FNAME = "supp_"
LDOS_FNAME = "ldos_"
PSI_FNAME = "psi_"
EXT = ".dat"


class Cell(object):
	"""Simulation cell which holds Atom objects in a 3D mesh.

	All lengths measured in Bohr radii (a0).
	All energies measured in Hartrees (Ha).
	
	Attributes:
		name (string): Name of simulation; used for plot titles
			fermiLevel (float): Fermi Level of simulation
			xLength (float): Length of cell along x
			yLength (float): Length of cell along y
			zLength (float): Length of cell along z
			gridSpacing (float): Resolution of mesh points
			x_mesh (3D mgrid): Mesh of x values
			y_mesh (3D mgrid): Mesh of y values
			z_mesh (3D mgrid): Mesh of z values
			atoms (int : Atom): Atom objects of simulation indexed by atom number
			basis_points (SmartDict): Basis function values, indexed by [atomKey][x][y][z][l][zeta][m]
			bands (Vector : [float]): Energies of bands indexed by k-point vector
	"""

	def __init__(self, name, fermiLevel, xLength, yLength, zLength, grid_spacing=0.5):
		"""Constructs 3D cell with given dimensional.

		All lengths measured in Bohr radii (a0);
		All energies measured in Hartrees (Ha).

		Args:
			name (string): Name of simulation; used for plot titles
			fermiLevel (float): Fermi Level of simulation
			xLength (float): Length of cell along x
			yLength (float): Length of cell along y
			zLength (float): Length of cell along z
			grid_spacing (float, opt.): Resolution of mesh points
		"""

		self.name = name
		self.fermiLevel = fermiLevel
		self.grid_spacing = grid_spacing

		vector_x = int(xLength/grid_spacing)*grid_spacing
		vector_y = int(yLength/grid_spacing)*grid_spacing
		vector_z = int(zLength/grid_spacing)*grid_spacing
		self.vector = Vector(vector_x, vector_y, vector_z)

		# Form Cartesian meshes
		self.x_mesh, self.y_mesh, self.z_mesh = np.mgrid[0: xLength: grid_spacing, 0: yLength: grid_spacing, 0: zLength: grid_spacing]

		# Initialise atoms and bands
		self.atoms = {}
		self.bands = {}
		self.support_grid = None

	def has_band(self, K, E):
		"""Check if cell stores specified band.

		Args:
			K (Vector): 3D Cartesian k-space vector
			E (float): Band energy
		"""
		output = False
		if K in self.bands:
			if E in self.bands[K]:
				output = True
		return output

	# def constrain_relative_vector(self, vector):
	# 	"""Return a vector that is constrained within simulation cell"""
	# 	new_vector = deepcopy(vector)
	#
	# 	# Check if vector components are greater than half of cell sides
	# 	# If greater, add or subtract cell length
	#
	# 	if vector.x > self.vector.x/2:
	# 		new_vector -= self.vector.project_x()
	# 	elif vector.x <= -self.vector.x/2:
	# 		new_vector += self.vector.project_x()
	#
	# 	if vector.y > self.vector.y/2:
	# 		new_vector -= self.vector.project_y()
	# 	elif vector.y <= -self.vector.y/2:
	# 		new_vector += self.vector.project_y()
	#
	# 	if vector.z > self.vector.z/2:
	# 		new_vector -= self.vector.project_z()
	# 	elif vector.z <= -self.vector.z/2:
	# 		new_vector += self.vector.project_z()
	#
	# 	# If vector is unchanged return
	# 	# If vector has changed, constrain
	# 	if new_vector == vector:
	# 		return new_vector
	# 	else:
	# 		return self.constrain_relative_vector(new_vector)
	#
	# def constrain_vector_to_cell(self, vector):
	# 	"""Return a vector that is constrained within simulation cell"""
	# 	new_vector = deepcopy(vector)
	#
	# 	# Check if vector components are greater than half of cell sides
	# 	# If greater, add or subtract cell length
	#
	# 	if vector.x >= self.vector.x:
	# 		new_vector -= self.vector.project_x()
	# 	elif vector.x < 0:
	# 		new_vector += self.vector.project_x()
	#
	# 	if vector.y >= self.vector.y:
	# 		new_vector -= self.vector.project_y()
	# 	elif vector.y < 0:
	# 		new_vector += self.vector.project_y()
	#
	# 	if vector.z >= self.vector.z:
	# 		new_vector -= self.vector.project_z()
	# 	elif vector.z < 0:
	# 		new_vector += self.vector.project_z()
	#
	# 	# If vector is unchanged return
	# 	# If vector has changed, constrain
	# 	if new_vector == vector:
	# 		return new_vector
	# 	else:
	# 		return self.constrain_vector_to_cell(new_vector)

	def add_atom(self, atom, atomKey):
		"""Add atom to self.atoms, indexed by atomKey

		Args:
			atom (Atom): Atom object
			atomKey (int): Atom number, as given in Conquest_out
		"""
		# Add to dict
		self.atoms[atomKey] = atom
		# Loop over atoms k-points
		for K in atom.bands:
			# If cell does not have k-point, create empty band energy list
			if K not in self.bands:
				self.bands[K] = []
			# Add band energies to k-point
			for E in atom.bands[K]:
				if E not in self.bands[K]:
					self.bands[K].append(E)
			# Sort energy list
			self.bands[K] = sorted(self.bands[K])

	def fermi_dirac(self, energy, temperature):
		"""Calculate Fermi-Dirac distribution value.

		Args:
			energy (float): Energy in eV
			temperature (float): Absolute temperature in K

		Returns:
			float: Occupation value
		"""
		f = 1.0 / (np.exp((energy - self.fermiLevel) / (BOLTZMANN * temperature)) + 1)
		return f

	def get_nearest_mesh_value(self, x):
		"""Return nearest mesh point to x.

		Works for any direction: x, y, or z.
		"""
		# Get quotient and remainder wrt grid spacing
		div_x, mod_x = divmod(x, self.grid_spacing)
		# Check if x should be rounded up or down
		if mod_x >= self.grid_spacing / 2:
			div_x += 1
		# Get new point
		new_x = div_x * self.grid_spacing
		return new_x

	def get_nearest_mesh_vector(self, position):
		"""Return vector to nearest mesh point to position vector"""
		# Get nearest mesh point for each component
		x = self.get_nearest_mesh_value(position.x)
		y = self.get_nearest_mesh_value(position.y)
		z = self.get_nearest_mesh_value(position.z)
		return Vector(x, y, z)

	def constrain_relative_vector(self, vector):
		"""Return a vector that is constrained within simulation cell"""
		x, y, z = vector.components()

		# Check if vector components are greater than half of cell sides
		# If greater, add or subtract cell length

		while x > self.vector.x / 2:
			x -= self.vector.x
		while x <= -self.vector.x / 2:
			x += self.vector.x

		while y > self.vector.y / 2:
			y -= self.vector.y
		while y <= -self.vector.y / 2:
			y += self.vector.y

		while z > self.vector.z / 2:
			z -= self.vector.z
		while z <= -self.vector.z / 2:
			z += self.vector.z

		return Vector(x, y, z)

	def constrain_vector_to_cell(self, vector):
		"""Return a vector that is constrained within simulation cell"""
		x, y, z = vector.components()

		# Check if vector components are greater than half of cell sides
		# If greater, add or subtract cell length

		while x >= self.vector.x:
			x -= self.vector.x
		while x < 0:
			x += self.vector.x

		while y >= self.vector.y:
			y -= self.vector.y
		while y < 0:
			y += self.vector.y

		while z >= self.vector.z:
			z -= self.vector.z
		while z < 0:
			z += self.vector.z

		return Vector(x, y, z)

	def calculate_support_grid(self, debug=False):
		"""Evaluate support function for each PAO on mesh.

		Args:
			debug (bool, opt.): If true, print extra information during runtime

		Returns:
			3D np array of SmartDict: Support function mesh, indexed by [x, y, z][atom_key][l][zeta][m]
		"""
		if debug:
			print "Calculating support grid"
		# Initialise support grid
		support_grid = np.empty_like(self.x_mesh, dtype=SmartDict)

		# Iterate over all atoms
		for atom_key in self.atoms:
			atom = self.atoms[atom_key]

			# Get atom cutoff radius
			cut = atom.get_max_cutoff()

			# Get nearest mesh point to atom position
			atom_pos_on_mesh = self.get_nearest_mesh_vector(atom.atom_pos)

			# Get mesh points of maximum range of atoms orbitals in each direction
			x_lower_lim = self.get_nearest_mesh_value(atom.atom_pos.x - cut)
			x_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.x + cut) + self.grid_spacing
			y_lower_lim = self.get_nearest_mesh_value(atom.atom_pos.y - cut)
			y_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.y + cut) + self.grid_spacing
			z_lower_lim = self.get_nearest_mesh_value(atom.atom_pos.z - cut)
			z_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.z + cut) + self.grid_spacing

			# Get array of mesh points within cutoff
			x_points = np.arange(x_lower_lim, x_upper_lim, self.grid_spacing)
			y_points = np.arange(y_lower_lim, y_upper_lim, self.grid_spacing)
			z_points = np.arange(z_lower_lim, z_upper_lim, self.grid_spacing)

			# Iterate over mesh points within cutoff
			for x in x_points:
				for y in y_points:
					for z in z_points:
						r = Vector(x, y, z)
						# Constrain vector using periodic boundary conditions
						constrained = self.constrain_vector_to_cell(r)

						# Get indices of cell mesh corresponding to this point
						i = np.where(self.x_mesh == constrained.x)[0][0]
						j = np.where(self.y_mesh == constrained.y)[1][0]
						k = np.where(self.z_mesh == constrained.z)[2][0]

						relative_position = self.constrain_relative_vector(r - atom_pos_on_mesh)
						# Iterate over orbitals
						for l in atom.radials:
							for zeta in atom.radials[l]:
								# Get radial part of wavefunction
								R = atom.get_radial_value_relative(l, zeta, relative_position)
								# If R == 0, do not store
								if R != 0.0:
									for m in range(-l, l + 1):
										# Get spherical harmonic
										Y = sph(l, m, relative_position)
										# Initialise support grid entry
										if not support_grid[i, j, k]:
											support_grid[i, j, k] = SmartDict()
										# Store support value
										support_grid[i, j, k][atom_key][l][zeta][m] = R * Y
			if debug:
				print "Calculated support grid for atom " + atom.ion_name + " " + str(atom_key)
		return support_grid

	def support_filename(self):
		return MESH_FOLDER+SUPPORT_FNAME+self.name+"_"+str(self.grid_spacing)+EXT

	def has_support_file(self):
		output = False
		filename = self.support_filename()
		if os.path.isfile(filename):
			support_file = open(filename, 'r')
			try:
				line = support_file.next()
				spacing_string = line.split()[0]
				if float(spacing_string) == self.grid_spacing:
					output = True
			except StopIteration:
				pass
			except ValueError:
				pass
			support_file.close()
		return output

	def write_support_grid(self, debug=False):
		"""Write support function mesh to file.

		Args:
			debug (bool, opt.): Print extra information during runtime
		"""
		filename = self.support_filename()
		support_file = safe_open(filename, 'w')
		if debug:
			print "Writing support grid to "+filename

		# Iterate over mesh points
		for i in range(self.x_mesh.shape[0]):
			for j in range(self.y_mesh.shape[1]):
				for k in range(self.z_mesh.shape[2]):
					# If support function values exist at mesh point
					if self.support_grid[i, j, k]:
						support_file.write(str(i)+" "+str(j)+" "+str(k)+"\n")
						# Iterate over atoms
						for atom_key in self.support_grid[i, j, k]:
							# Write atom index
							support_file.write(str(atom_key)+"\n")
							# Iterate over orbitals
							for l in self.support_grid[i, j, k][atom_key]:
								for zeta in self.support_grid[i, j, k][atom_key][l]:
									for m in self.support_grid[i, j, k][atom_key][l][zeta]:
										# Write orbital data
										line = (str(l)+" "+str(zeta)+" "+str(m)+" "
												+str(self.support_grid[i, j, k][atom_key][l][zeta][m]))
										support_file.write(line+"\n")
		support_file.close()

	def read_support_grid(self, debug=False):
		"""Read support function mesh from file"""
		filename = self.support_filename()
		support_file = open(filename, 'r')
		support_grid = np.empty_like(self.x_mesh, dtype=SmartDict)

		if debug:
			print "Reading support grid from "+filename

		# Iterate over file lines
		end_of_file = False
		line = support_file.next()
		line_split = line.split()
		while not end_of_file:
			try:
				# Get mesh indices
				i = int(line_split[0])
				j = int(line_split[1])
				k = int(line_split[2])
				# Read atom data
				reading_atoms = True
				line = support_file.next()
				line_split = line.split()
				while reading_atoms:
					if len(line_split) == 1:
						# Get atom index
						atom_key = int(line)
						# Read orbital data
						reading_orbitals = True
						while reading_orbitals:
							line = support_file.next()
							line_split = line.split()
							if len(line_split) != 4:
								reading_orbitals = False
							elif len(line_split) == 1:
								reading_atoms = True
							else:
								l = int(line_split[0])
								zeta = int(line_split[1])
								m = int(line_split[2])
								value = float(line_split[3])
								if not support_grid[i, j, k]:
									support_grid[i, j, k] = SmartDict()
								support_grid[i, j, k][atom_key][l][zeta][m] = value
					else:
						reading_atoms = False
			except StopIteration:
				end_of_file = True
		if debug:
			print "Support grid successfully read"
		return support_grid

	def get_support_grid(self, recalculate=False, write=True, debug=False):
		"""Get support function mesh.

		Args:
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file
			debug (bool, opt.): Print extra information during runtime

		Returns:
			3D np.array of SmartDict: Support function mesh, indexed by [x, y, z][atomKey][l][zeta][m]
		"""
		if self.support_grid is None:
			if not recalculate and os.path.isfile(self.support_filename()):
				# Read support grid from file
				self.support_grid = self.read_support_grid(debug=debug)
			else:
				# Recalculate support grid
				self.support_grid = self.calculate_support_grid(debug=debug)
				# Write to file
				if write:
					self.write_support_grid(debug=debug)
		return self.support_grid

	def calculate_psi_grid(self, K, E, recalculate=False, write=True, debug=False):
		if debug:
			print "Building Wavefunction for "+str(K)+", "+str(E)
		psi_grid = np.zeros_like(self.x_mesh, dtype=complex)
		support_grid = self.get_support_grid(recalculate=recalculate, write=write, debug=debug)
		for atom_key in self.atoms:
			atom = self.atoms[atom_key]
			for l in atom.bands[K][E]:
				for zeta in atom.bands[K][E][l]:
					for m in atom.bands[K][E][l][zeta]:
						coefficient = atom.get_coefficient(K, E, l, zeta, m)
						for i in range(self.x_mesh.shape[0]):
							for j in range(self.y_mesh.shape[1]):
								for k in range(self.z_mesh.shape[2]):
									if support_grid[i, j, k]:
										if atom_key in support_grid[i, j, k]:
											psi_grid[i, j, k] += coefficient*support_grid[i, j, k][atom_key][l][zeta][m]
		return psi_grid

	def psi_filename(self, K, E):
		return (MESH_FOLDER+PSI_FNAME+self.name+"_"+str(self.grid_spacing)+"_"
				+str(K.x)+"_"+str(K.y)+"_"+str(K.z)+"_"+str(E)+EXT)

	def write_psi_grid(self, psi_grid, K, E):
		filename = self.psi_filename(K, E)
		psi_file = safe_open(filename, "w")
		for i in range(self.x_mesh.shape[0]):
			for j in range(self.y_mesh.shape[1]):
				for k in range(self.z_mesh.shape[2]):
					if psi_grid[i, j, k] != 0:
						psi = psi_grid[i, j, k]
						psi_file.write(str(i)+" "+str(j)+" "+str(k)+" "+str(psi.real)+" "+str(psi.imag)+"\n")
		psi_file.close()

	def write_psi_range(self, min_E, max_E, recalculate=False, debug=False):
		for K in self.bands:
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					psi_grid = self.calculate_psi_grid(K, E, recalculate=recalculate, debug=debug)
					self.write_psi_grid(psi_grid, K, E)

	def read_psi_grid(self, K, E):
		filename = self.psi_filename(K, E)
		psi_file = open(filename, "r")
		psi_grid = np.zeros_like(self.x_mesh, dtype=complex)

		for line in psi_file:
			line_split = line.split()
			i = int(line_split[0])
			j = int(line_split[1])
			k = int(line_split[2])
			psi_real = float(line_split[3])
			psi_imag = float(line_split[4])
			psi = complex(psi_real, psi_imag)
			psi_grid[i, j, k] = psi
		return psi_grid

	def get_psi_grid(self, K, E, recalculate=False, write=True, debug=False):
		"""Get mesh of complex wavefunction values.

		Args:
			K (KVector): 3D Cartesian k-vector
			E (float): Band energy
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file
			debug (bool, opt.): Print extra information during runtime
		"""
		if not recalculate and os.path.isfile(self.psi_filename(K, E)):
			# Read data from file
			psi_grid = self.read_psi_grid(K, E)
		else:
			psi_grid = self.calculate_psi_grid(K, E, recalculate=recalculate, write=write, debug=debug)
			if write:
				self.write_psi_grid(psi_grid, K, E)
		return psi_grid

	def calculate_ldos_grid(self, min_E, max_E, T, recalculate=False, write=True, debug=False):
		"""Calculate LDoS mesh.

		Args:
			min_E: Minimum energy
			max_E: Maximum energy
			T: Absolute temperature in Kelvin
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write calculated meshes to file
			debug (bool, opt.): Print extra information during runtime

		Returns:
			3D np.array: LDoS mesh
		"""

		if debug:
			print "Calculating LDoS grid"
		ldos_grid = np.zeros_like(self.x_mesh, dtype=float)

		for K in self.bands:
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					psi_grid = self.get_psi_grid(K, E, recalculate=recalculate, write=write, debug=debug)
					fd = self.fermi_dirac(E, T)
					for i in range(self.x_mesh.shape[0]):
						for j in range(self.y_mesh.shape[1]):
							for k in range(self.z_mesh.shape[2]):
								psi = psi_grid[i, j, k]
								ldos_grid[i, j, k] += K.weight*fd*(abs(psi))**2
				if debug:
					print "Completed LDoS for ", K, E
		return ldos_grid

	def ldos_filename(self, min_E, max_E, T):
		return MESH_FOLDER+LDOS_FNAME+self.name+"_"+str(self.grid_spacing)+"_"+str(min_E)+"_"+str(max_E)+"_"+str(T)+EXT

	def write_ldos_grid(self, ldos_grid, min_E, max_E, T, debug=False):
		"""Write LDoS mesh to file.

		Args:
			min_E: Minimum energy
			max_E: Maximum energy
			T: Absolute temperature in K
			debug (bool, opt.): Print extra information during runtime
		"""
		filename = self.ldos_filename(min_E, max_E, T)
		# Get LDoS mesh
		ldos_file = safe_open(filename, "w")
		if debug:
			print "Writing LDoS grid to "+filename
		# Iterate over mesh points
		for i in range(self.x_mesh.shape[0]):
			for j in range(self.y_mesh.shape[1]):
				for k in range(self.z_mesh.shape[2]):
					# If LDoS is non-zero at mesh point, write data to file
					if ldos_grid[i, j, k]:
						ldos_file.write(str(i)+" "+str(j)+" "+str(k)+" "+str(ldos_grid[i, j, k])+"\n")
		ldos_file.close()

	def read_ldos_grid(self, min_E, max_E, T, debug=False):
		"""Read LDoS mesh from file"""
		filename = self.ldos_filename(min_E, max_E, T)
		ldos_file = open(filename, 'r')
		ldos_grid = np.zeros_like(self.x_mesh)

		if debug:
			print "Reading LDoS grid from "+filename

		for line in ldos_file:
			line_split = line.split()
			# Get mesh indices
			i = int(line_split[0])
			j = int(line_split[1])
			k = int(line_split[2])
			# Get LDoS value
			value = float(line_split[3])
			ldos_grid[i, j, k] = value

		if debug:
			print "LDoS grid successfully read"
		return ldos_grid

	def get_ldos_grid(self, min_E, max_E, T, recalculate=False, write=True, debug=False):
		"""Get LDoS mesh.

		Args:
			min_E: Minimum energy
			max_E: Maximum energy
			T: Absolute temperature in K
			recalculate (bool, opt.): Force recalculation of meshes, even if already stored
			write (bool, opt.): Write calculated grids to file
			debug (bool, opt.): Print extra information during runtime

		Returns:
			3D np.array: LDoS mesh
		"""
		# Read ldos grid from file if not stored by cell
		print self.ldos_filename(min_E, max_E, T)
		if not recalculate and os.path.isfile(self.ldos_filename(min_E, max_E, T)):
			ldos_grid = self.read_ldos_grid(min_E, max_E, T, debug=debug)
		else:
			# Calculate LDoS on mesh
			ldos_grid = self.calculate_ldos_grid(min_E, max_E, T, recalculate=recalculate, write=write, debug=debug)
			if write:
				self.write_ldos_grid(ldos_grid, min_E, max_E, T, debug=debug)
		return ldos_grid

	def kappa_squared(self, tip_work_func, tip_energy):
		return 2*ELECTRON_MASS/(H_BAR**2)*(tip_work_func - tip_energy)

	def greens_function(self, position, tip_position, tip_work_func, tip_energy):
		distance = abs(position - tip_position)
		kappa2 = self.kappa_squared(tip_work_func, tip_energy)
		return np.exp(- kappa2 * distance) / (4*np.pi*distance)

	def filtered_greens_function(self, position, tip_position, alpha, tip_work_func, tip_energy):
		distance = abs(position - tip_position)
		kappa2 = self.kappa_squared(tip_work_func, tip_energy)
		kappa = np.sqrt(kappa2)
		u = np.sqrt(2*np.pi)/(4*distance)*np.exp(alpha**2*kappa2)
		v = np.exp(np.sqrt(kappa)*distance)*math.erf(alpha*kappa + distance/(2*alpha))
		w = np.exp(- kappa * distance)*math.erf(alpha*kappa - distance/(2*alpha))
		x = 2*np.sinh(kappa * distance)
		return u*(v - w - x)

	def isosurface_condition(self, ldos, reference):
		return np.log(ldos / reference)

	def broadened_volume_integrand(self, mesh, width, fraction):
		masked_mesh = np.ma.masked_array(mesh, mask=np.where(mesh == 0, 1, 0))
		log_mesh = self.isosurface_condition(masked_mesh, fraction*np.max(mesh))
		broad_mesh = np.where(log_mesh < width, 15/(16*width)*(1-(log_mesh/width)**2)**2, 0)
		return broad_mesh
