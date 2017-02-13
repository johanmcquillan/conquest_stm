
import os
import sys
import warnings
import numpy as np
from sph import sph
from smart_dict import SmartDict
from vector import Vector, KVector
from safe_open import safe_open


class Cell(object):
	"""Simulation cell which holds Atom objects in a 3D mesh.

	All lengths measured in Bohr radii (a0).
	All energies measured in electron volts (eV).
	
	Attributes:
		name (string): Name of simulation; used for plot titles
		fermi_level (float): Fermi Level of simulation
		vector (Vector): Vector from origin to far corner of cell; components are maximum lengths of cell
		grid_spacing (float): Resolution of mesh points
		real_mesh (np.ndarray): Real space mesh
		atoms ({int : Atom}): Atom objects of simulation indexed by atom number
		bands ({Vector : [float]}): Energies of bands indexed by k-point vector
		support_mesh (array(SmartDict)): Mesh of support function values, indexed by [x, y, z][atom_key][l][zeta][m]
	"""

	BOLTZMANN = 8.6173303E-5  # Boltzmann's Constant in eV/K
	ELECTRON_MASS = 9.10938E-31  # Electron Mass in kg
	H_BAR = 4.135667662E-15  # Reduced Planck's Constant in eV.s
	DELTA_S_FACTOR = 1.5

	MESH_FOLDER = "meshes/"
	SUPPORT_FNAME = "supp_"
	LDOS_FNAME = "ldos_"
	PSI_FNAME = "psi_"
	PROP_PSI_FNAME = "prop_"
	CURRENT_FNAME = "current_"
	CURRENT_PLANE_FNAME = "cur_plane_"
	EXT = ".dat"

	PRINT_RELATIVE_TO_EF = True
	PROG_BAR_INTERVALS = 20
	PROG_BAR_CHARACTER = ">"

	def __init__(self, name, fermi_level, x_length, y_length, z_length, grid_spacing=0.5, group_size=150):
		"""Constructs 3D cell with given dimensional.

		All lengths measured in Bohr radii (a0);
		All energies measured in Hartrees (Ha).

		Args:
			name (string): Name of simulation; used for plot titles
			fermi_level (float): Fermi Level of simulation
			x_length (float): Length of cell along x
			y_length (float): Length of cell along y
			z_length (float): Length of cell along z
			grid_spacing (float, opt.): Resolution of mesh points
		"""

		self.name = name
		self.fermi_level = fermi_level
		self.grid_spacing = grid_spacing
		self.atom_group_size = group_size

		vector_x = int(x_length / grid_spacing) * grid_spacing
		vector_y = int(y_length / grid_spacing) * grid_spacing
		vector_z = int(z_length / grid_spacing) * grid_spacing
		self.vector = Vector(vector_x, vector_y, vector_z)

		# Form Cartesian meshes
		self.real_mesh = np.transpose(np.mgrid[0:x_length:grid_spacing,
		                                       0: y_length: grid_spacing, 0:
		                                       z_length: grid_spacing],
		                              (1, 2, 3, 0))
		self.mesh_points = self.real_mesh.shape[0] * self.real_mesh.shape[1] * self.real_mesh.shape[2]

		# Initialise atoms and bands
		self.atoms = {}
		self.bands = {}
		self.support_mesh = None
		self.current_group = -1
		self.default_delta_s = self.delta_s()

		self.psi_vec = np.vectorize(self.calculate_psi_grid_vec)

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

	def add_atom(self, atom, atomKey):
		"""Add atom to self.atoms, indexed by atom_key

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
		"""Calculate Fermi-Dirac occupation factor.

		Args:
			energy (float): Energy in eV
			temperature (float): Absolute temperature in K

		Returns:
			float: Occupation factor, between 0 and 1
		"""
		with warnings.catch_warnings():
			# Suppress RuntimeWarning from overflow and underflow in np.exp
			warnings.simplefilter('ignore', RuntimeWarning)
			f = 1.0 / (np.exp((energy - self.fermi_level) / (self.BOLTZMANN * temperature)) + 1)
		return f

	def get_nearest_mesh_value(self, x, indices=False, points=None):
		"""Return nearest mesh point to x. Not constrained within simulation cell.

		Works for any direction: x, y, or z.
		"""
		# Get quotient and remainder wrt grid spacing
		div_x, mod_x = divmod(x, self.grid_spacing)
		# Check if x should be rounded up or down
		if mod_x >= self.grid_spacing / 2:
			div_x += 1
		# Get new point
		new_x = div_x * self.grid_spacing
		if indices and points is not None:
			i = div_x - 1
			while i < 0:
				i += points
			while i >= points:
				i -= points
			return new_x, int(i)
		else:
			return new_x

	def get_nearest_mesh_vector(self, position, indices=False):
		"""Return vector to nearest mesh point to position vector"""
		# Get nearest mesh point for each component
		if indices:
			x, i = self.get_nearest_mesh_value(position.x, indices=True)
			y, j = self.get_nearest_mesh_value(position.y, indices=True)
			z, k = self.get_nearest_mesh_value(position.z, indices=True)
			ijk = (i, j, k)
			return Vector(x, y, z), ijk
		else:
			x = self.get_nearest_mesh_value(position.x, indices=False)
			y = self.get_nearest_mesh_value(position.y, indices=False)
			z = self.get_nearest_mesh_value(position.z, indices=False)
			return Vector(x, y, z)

	def constrain_relative_vector(self, vector):
		"""Return a vector that is constrained within simulation cell"""
		x, y, z = vector.components

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
		x, y, z = vector.components

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

	def calculate_support_group(self, group, debug=False):
		"""Evaluate support function for each PAO on mesh.

		Args:
			vectorised (bool, opt.): If true, use NumPy vectorisation
			debug (bool, opt.): If true, print extra information during runtime

		Returns:
			array(SmartDict): Mesh of support function values, indexed by [x, y, z][atom_key][l][zeta][m]
		"""
		if debug:
			print "Calculating support grid"
		# Initialise support grid
		support_grid = np.empty(self.real_mesh.shape[:3], dtype=SmartDict)

		# Iterate over all atoms
		lower_bound = group * self.atom_group_size
		upper_bound = (group + 1) * self.atom_group_size
		atom_list = sorted([a for a in self.atoms.keys() if lower_bound <= a < upper_bound])
		for atom_key in atom_list:
			atom = self.atoms[atom_key]

			debug_str = "  Support grid for atom {:3} - {!s}: ".format(atom_key, atom.ion_name)
			if debug:
				sys.stdout.write(debug_str)
				sys.stdout.flush()

			# Get atom cutoff radius
			cut = atom.get_max_cutoff()

			# Get nearest mesh point to atom position
			atom_pos_on_mesh = self.get_nearest_mesh_vector(atom.atom_pos)

			# Get mesh points of maximum range of atoms orbitals in each direction
			x_lower_lim, i_start = self.get_nearest_mesh_value(atom.atom_pos.x - cut, indices=True, points=self.real_mesh.shape[0])
			x_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.x + cut) + self.grid_spacing
			y_lower_lim, j_start = self.get_nearest_mesh_value(atom.atom_pos.y - cut, indices=True, points=self.real_mesh.shape[1])
			y_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.y + cut) + self.grid_spacing
			z_lower_lim, k_start = self.get_nearest_mesh_value(atom.atom_pos.z - cut, indices=True, points=self.real_mesh.shape[2])
			z_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.z + cut) + self.grid_spacing

			# Get array of mesh points within cutoff
			local_mesh = np.transpose(np.mgrid[x_lower_lim:x_upper_lim:self.grid_spacing, y_lower_lim:y_upper_lim:self.grid_spacing, z_lower_lim:z_upper_lim:self.grid_spacing], (1, 2, 3, 0))
			lm_shape = local_mesh.shape[:3]
			points_done = 0
			bars_done = 0
			total_points = local_mesh.shape[0]*local_mesh.shape[1]*local_mesh.shape[2]

			rolled_mesh = np.roll(support_grid, -i_start, 0)
			rolled_mesh = np.roll(rolled_mesh, -j_start, 1)
			rolled_mesh = np.roll(rolled_mesh, -k_start, 2)
			partial_mesh = rolled_mesh[0:lm_shape[0], 0:lm_shape[1], 0:lm_shape[2]]

			for local_ijk in np.ndindex(lm_shape):
				position = local_mesh[local_ijk]
				r = Vector(*position)
				relative_position = self.constrain_relative_vector(r - atom_pos_on_mesh)

				partial_ijk_list = list(local_ijk)
				for i in range(len(partial_ijk_list)):
					if partial_ijk_list[i] >= self.real_mesh.shape[i]:
						partial_ijk_list[i] -= self.real_mesh.shape[i]
				partial_ijk = tuple(partial_ijk_list)

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
								if partial_mesh[partial_ijk] is None:
									partial_mesh[partial_ijk] = SmartDict()

								if m not in partial_mesh[partial_ijk][atom_key][l][zeta]:
									partial_mesh[partial_ijk][atom_key][l][zeta][m] = 0
								# Store support value
								partial_mesh[partial_ijk][atom_key][l][zeta][m] += R * Y
				points_done += 1

				prog = float(points_done) / total_points
				if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
					percent = prog * 100
					sys.stdout.write('\r')
					sys.stdout.write(debug_str)
					sys.stdout.write(
						" [{:<{}}]".format(self.PROG_BAR_CHARACTER * bars_done, self.PROG_BAR_INTERVALS))
					sys.stdout.write(" {:3.0f}%".format(percent))
					sys.stdout.flush()
					bars_done += 1

			rolled_mesh[0:lm_shape[0], 0:lm_shape[1], 0:lm_shape[2]] = partial_mesh
			rolled_mesh = np.roll(rolled_mesh, i_start, 0)
			rolled_mesh = np.roll(rolled_mesh, j_start, 1)
			support_grid = np.roll(rolled_mesh, k_start, 2)

			if debug:
				sys.stdout.write("\n")
				sys.stdout.flush()

		return support_grid

	def support_group_filename(self, group):
		"""Return standardised filename for relevant support function file"""
		return self.MESH_FOLDER+self.SUPPORT_FNAME+self.name+"_"+str(self.grid_spacing)+"_"+str(self.atom_group_size) + "_" + str(group) + self.EXT

	def write_support_group(self, group, support_mesh, debug=False):
		"""Write support function mesh to file"""
		filename = self.support_group_filename(group)
		support_file = safe_open(filename, 'w')

		if debug:
			sys.stdout.write("Writing support group to "+filename+": ")
			sys.stdout.flush()
		points_done = 0
		bars_done = 0

		# Iterate over mesh points
		for ijk in np.ndindex(self.real_mesh.shape[:3]):
			i, j, k = ijk

			# If support function values exist at mesh point
			if support_mesh[ijk]:
				support_file.write(str(i) + " " + str(j) + " " + str(k) + "\n")

				# Iterate over atoms
				for atom_key in support_mesh[ijk]:
					# Write atom index
					support_file.write(str(atom_key) + "\n")
					# Iterate over orbitals

					for l in support_mesh[ijk][atom_key]:
						for zeta in support_mesh[ijk][atom_key][l]:
							for m in support_mesh[ijk][atom_key][l][zeta]:
								# Write orbital data
								line = (str(l) + " " + str(zeta) + " " + str(m) + " "
								        + str(support_mesh[ijk][atom_key][l][zeta][m]))
								support_file.write(line + "\n")

			points_done += 1
			if debug and float(points_done) / self.mesh_points * self.PROG_BAR_INTERVALS > bars_done:
				sys.stdout.write(" [{:<{}}]".format(self.PROG_BAR_CHARACTER * bars_done, self.PROG_BAR_INTERVALS))
				sys.stdout.flush()
				bars_done += 1

		if debug:
			sys.stdout.write("\n")
			sys.stdout.flush()
		support_file.close()

	def read_support_group(self, group, debug=False):
		"""Read support function mesh from file"""
		filename = self.support_group_filename(group)
		support_file = open(filename, 'r')
		support_mesh = np.empty(self.real_mesh.shape[:3], dtype=SmartDict)

		if debug:
			print "Reading support group from "+filename

		# Iterate over file lines
		end_of_file = False
		line = support_file.next()
		line_split = line.split()
		while not end_of_file:
			try:
				# Get mesh indices
				i, j, k = [int(a) for a in line_split[:3]]
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
								if support_mesh[i, j, k] is None:
									support_mesh[i, j, k] = SmartDict()
								support_mesh[i, j, k][atom_key][l][zeta][m] = value
					else:
						reading_atoms = False
			except StopIteration:
				end_of_file = True
		if debug:
			print "Support grid successfully read"
		return support_mesh

	def get_support_group(self, group, recalculate=False, write=True, vectorised=True, debug=False):
		"""Get support function mesh.

		Args:
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file
			vectorised (bool, opt.): If true, use NumPy vectorisation
			debug (bool, opt.): Print extra information during runtime

		Returns:
			array(SmartDict): Mesh of support function values, indexed by [x, y, z][atom_key][l][zeta][m]
		"""
		if group == self.current_group and self.support_mesh is not None:
			return self.support_mesh
		else:
			if not recalculate and os.path.isfile(self.support_group_filename(group)):
				# Read support grid from file
				support_mesh = self.read_support_group(group, debug=debug)
			else:
				# Recalculate support grid
				support_mesh = self.calculate_support_group(group, debug=debug)
				# Write to file
				if write:
					self.write_support_group(group, support_mesh, debug=debug)
			self.current_group = group
			self.support_mesh = support_mesh
			return support_mesh

	def calculate_psi_grid_vec(self, support_dict, atom_key, l, zeta, m, coefficient):
		# Evaluate wavefunction contribution
		if support_dict is not None and atom_key in support_dict and l in support_dict[atom_key] and zeta in support_dict[atom_key][l] and m in support_dict[atom_key][l][zeta]:
			psi = coefficient * support_dict[atom_key][l][zeta][m]
		else:
			psi = complex(0, 0)
		return psi

	def calculate_psi_grid(self, K, E, recalculate=False, write=True, vectorised=True, debug=False):
		"""Evaluate wavefunction on mesh.

		Args:
			K (KVector): K-point
			E (float): Band energy
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file
			vectorised (bool, opt.): If true, use NumPy vectorisation
			debug (bool, opt.): Print extra information during runtime

		Returns:
			array(complex): Mesh of complex wavefunction values
		"""
		atoms_done = 0
		total_atoms = len(self.atoms)
		bars_done = 0

		# Initialise mesh
		psi_grid = np.zeros_like(self.real_mesh[..., 0], dtype=complex)

		# Get basis functions
		# support_grid = self.get_support_grid(recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)

		if self.PRINT_RELATIVE_TO_EF:
			E_str = str(E - self.fermi_level) + " eV"
		else:
			E_str = str(E) + " eV"
		debug_str = "Calculating psi(r) at k = "+str(K)+", E = "+E_str+": "

		if debug:
			sys.stdout.write(debug_str)
			sys.stdout.flush()

		previous_group = 0
		support_grid = self.get_support_group(0, debug=debug)
		for atom_key in sorted(self.atoms.iterkeys()):
			atom = self.atoms[atom_key]
			group = atom_key / self.atom_group_size
			if group != previous_group:
				support_grid = self.get_support_group(group, debug=debug)
				previous_group = group
			# Iterate over orbitals
			for l in atom.bands[K][E]:
				for zeta in atom.bands[K][E][l]:
					for m in atom.bands[K][E][l][zeta]:
						# Evaluate wavefunction contribution over mesh
						coefficient = atom.get_coefficient(K, E, l, zeta, m)
						if vectorised:
							psi_grid += self.psi_vec(support_grid, atom_key, l, zeta, m, coefficient)
						else:
							for ijk in np.ndindex(self.real_mesh.shape[:3]):
								if support_grid[ijk]:
									if atom_key in support_grid[ijk]:
										psi_grid[ijk] += coefficient*support_grid[ijk][atom_key][l][zeta][m]
			# Print progress bar
			atoms_done += 1
			prog = float(atoms_done) / total_atoms
			if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
				percent = prog * 100
				sys.stdout.write('\r')
				sys.stdout.write(debug_str)
				sys.stdout.write(" [{:<{}}]".format(self.PROG_BAR_CHARACTER * bars_done, self.PROG_BAR_INTERVALS))
				sys.stdout.write(" {:3.0f}%".format(percent))
				sys.stdout.flush()
				bars_done += 1
		if debug:
			sys.stdout.write("\n")
			sys.stdout.flush()
		return psi_grid

	def psi_filename(self, K, E):
		"""Return standardised filename for relevant wavefunction file"""
		return (self.MESH_FOLDER+self.PSI_FNAME+self.name+"_"+str(self.grid_spacing)+"_"
				+str(K.x)+"_"+str(K.y)+"_"+str(K.z)+"_"+str(E)+self.EXT)

	def write_psi_grid(self, psi_grid, K, E):
		"""Write wavefunction function mesh to file"""
		filename = self.psi_filename(K, E)
		psi_file = safe_open(filename, "w")
		for ijk in np.ndindex(self.real_mesh.shape[:3]):
			i, j, k = ijk
			if psi_grid[ijk] != 0:
				psi = psi_grid[ijk]
				psi_file.write(str(i)+" "+str(j)+" "+str(k)+" "+str(psi.real)+" "+str(psi.imag)+"\n")
		psi_file.close()

	def read_psi_grid(self, K, E, debug=False, debug_file=False):
		"""Read wavefunction mesh from file"""
		filename = self.psi_filename(K, E)
		psi_file = open(filename, "r")
		psi_grid = np.zeros_like(self.real_mesh[..., 0], dtype=complex)

		if debug_file:
			sys.stdout.write("Reading psi(r) from {}\n".format(filename))
			sys.stdout.flush()
		elif debug:
			if self.PRINT_RELATIVE_TO_EF:
				E_str = str(E - self.fermi_level) + " eV"
			else:
				E_str = str(E) + " eV"
			sys.stdout.write("Reading psi(r) at k = {!s}, E = {}\n".format(K, E_str))

		for line in psi_file:
			line_split = line.split()
			i = int(line_split[0])
			j = int(line_split[1])
			k = int(line_split[2])
			real = float(line_split[3])
			imag = float(line_split[4])
			psi_grid[i, j, k] = complex(real, imag)
		return psi_grid

	def get_psi_grid(self, K, E, recalculate=False, write=True, vectorised=True, debug=False, debug_file=False):
		"""Get mesh of complex wavefunction values.

		Args:
			K (KVector): 3D Cartesian k-vector
			E (float): Band energy
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file
			vectorised (bool, opt.): If true, use NumPy vectorisation
			debug (bool, opt.): Print extra information during runtime
			debug_file (bool, opt.): If true, print name of files during debug rather than metadata

		Returns:
			array(complex): Mesh of complex wavefunction values
		"""
		if not recalculate and os.path.isfile(self.psi_filename(K, E)):
			# Read data from file
			psi_grid = self.read_psi_grid(K, E, debug=debug, debug_file=debug_file)
		else:
			psi_grid = self.calculate_psi_grid(K, E, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)
			if write:
				self.write_psi_grid(psi_grid, K, E)
		return psi_grid

	def calculate_ldos_grid(self, min_E, max_E, T, recalculate=False, write=True, vectorised=False, debug=False, debug_file=False):
		"""Calculate LDOS mesh.

		Args:
			min_E: Minimum energy
			max_E: Maximum energy
			T: Absolute temperature in Kelvin
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write calculated meshes to file
			vectorised (bool, opt.): If true, use NumPy vectorisation
			debug (bool, opt.): Print extra information during runtime
			debug_file (bool, opt.): If true, print name of files during debug rather than metadata

		Returns:
			3D np.array: LDoS mesh
		"""

		if debug:
			sys.stdout.write("Calculating local density of states grid\n")
			sys.stdout.flush()
		ldos_grid = np.zeros_like(self.real_mesh[..., 0], dtype=float)

		total_k_weight = 0
		for K in self.bands:
			total_k_weight += K.weight

		for K in self.bands:
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					psi_grid = self.get_psi_grid(K, E, recalculate=recalculate, write=write, debug=debug, debug_file=debug_file)
					fd = self.fermi_dirac(E, T)
					if E > self.fermi_level:
						fd = 1 - fd
					if vectorised:
						ldos_grid += (K.weight / total_k_weight) * fd * (abs(psi_grid))**2
					else:
						for ijk in np.ndindex(self.real_mesh.shape[:3]):
							ldos_grid[ijk] += (K.weight/total_k_weight)*fd*(abs(psi_grid[ijk]))**2
		return ldos_grid

	def ldos_filename(self, min_E, max_E, T):
		"""Return standardised filename for relevant LDOS file"""
		return self.MESH_FOLDER+self.LDOS_FNAME+self.name+"_"+str(self.grid_spacing)+"_"+str(min_E)+"_"+str(max_E)+"_"+str(T)+self.EXT

	def write_ldos_grid(self, ldos_grid, min_E, max_E, T, debug=False):
		"""Write LDoS mesh to file.

		Args:
			ldos_grid (array(float)): 3D array of LDOS values
			min_E (float): Minimum energy
			max_E (float): Maximum energy
			T (float): Absolute temperature in K
			debug (bool, opt.): Print extra information during runtime
		"""
		filename = self.ldos_filename(min_E, max_E, T)
		# Get LDOS mesh
		ldos_file = safe_open(filename, "w")
		if debug:
			sys.stdout.write("Writing LDOS grid to {}\n".format(filename))
			sys.stdout.flush()
		# Iterate over mesh points
		for i, j, k in np.ndindex(self.real_mesh.shape[:3]):
			# If LDOS is non-zero at mesh point, write data to file
			if ldos_grid[i, j, k]:
				ldos_file.write(str(i)+" "+str(j)+" "+str(k)+" "+str(ldos_grid[i, j, k])+"\n")
		ldos_file.close()

	def read_ldos_grid(self, min_E, max_E, T, debug=False):
		"""Read LDOS mesh from file"""
		filename = self.ldos_filename(min_E, max_E, T)
		ldos_file = open(filename, 'r')
		ldos_grid = np.zeros_like(self.real_mesh[..., 0])

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

	def get_ldos_grid(self, min_E, max_E, T, recalculate=False, write=True, vectorised=True, debug=False):
		"""Get LDOS mesh by calculating or reading from file.

		Args:
			min_E: Minimum energy
			max_E: Maximum energy
			T: Absolute temperature in K
			recalculate (bool, opt.): Force recalculation of meshes, even if already stored
			write (bool, opt.): Write calculated grids to file
			debug (bool, opt.): Print extra information during runtime

		Returns:
			array(float): 3D mesh of LDOS values
		"""
		# Read ldos grid from file if not stored by cell
		if not recalculate and os.path.isfile(self.ldos_filename(min_E, max_E, T)):
			ldos_grid = self.read_ldos_grid(min_E, max_E, T, debug=debug)
		else:
			# Calculate LDOS on mesh
			ldos_grid = self.calculate_ldos_grid(min_E, max_E, T, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)
			if write:
				self.write_ldos_grid(ldos_grid, min_E, max_E, T, debug=debug)
		return ldos_grid

	def get_vector_mesh(self, vectorised=True):
		"""Return 3D mesh of Vector objects pointing to corresponding element"""
		vector_mesh = np.empty(self.real_mesh.shape[:3], dtype=Vector)
		if vectorised:
			vector_mesh = np.vectorize(Vector)(self.real_mesh[..., 0], self.real_mesh[..., 1], self.real_mesh[..., 2])
		else:
			for ijk in np.ndindex(self.real_mesh.shape[:3]):
				x, y, z = self.real_mesh[ijk]
				vector_mesh[ijk] = Vector(x, y, z)
		return vector_mesh

	def periodic_gradient(self, mesh):
		"""Calculate gradient of mesh, enforcing periodic boundary conditions"""
		# Width of padding
		pad = 3
		# Shape of padded array
		padded_shape = (mesh.shape[0] + 2*pad, mesh.shape[1] + 2*pad, mesh.shape[2] + 2*pad)

		# Copy mesh into padded mesh
		padded_mesh = np.zeros(padded_shape, dtype=mesh.dtype)
		padded_mesh[pad:-pad, pad:-pad, pad:-pad] = mesh

		# Copy boundary regions into padding
		padded_mesh[:pad, pad:-pad, pad:-pad] = mesh[-pad:, :, :]
		padded_mesh[pad:-pad, :pad, pad:-pad] = mesh[:, -pad:, :]
		padded_mesh[pad:-pad, pad:-pad, :pad] = mesh[:, :, -pad:]
		padded_mesh[-pad:, pad:-pad, pad:-pad] = mesh[:pad, :, :]
		padded_mesh[pad:-pad, -pad:, pad:-pad] = mesh[:, :pad, :]
		padded_mesh[pad:-pad, pad:-pad, -pad:] = mesh[:, :, :pad]

		# Get gradient
		padded_gradient = np.transpose(np.array(np.gradient(padded_mesh, self.grid_spacing)), (1, 2, 3, 0))

		# Return unpadded gradient
		return padded_gradient[pad:-pad, pad:-pad, pad:-pad]

	def delta_s(self):
		min_delta_s = 2.0 * self.grid_spacing / self.H_BAR * np.sqrt(4.85 * 2.0 * C.ELECTRON_MASS)
		return self.DELTA_S_FACTOR * min_delta_s

	def kappa_squared(self, tip_work_func, tip_fermi_level):
		"""Calculate decay constant of tip wavefunction"""
		return 2*self.ELECTRON_MASS/(self.H_BAR**2)*(tip_work_func - tip_fermi_level)

	def greens_function(self, distance, tip_work_func, tip_energy):
		"""Evaluate Tersoff-Hamann Green's Function"""
		if distance == 0:
			return 0
		else:
			kappa2 = self.kappa_squared(tip_work_func, tip_energy)
			return np.exp(- kappa2 * distance) / (4*np.pi*distance)

	def broadened_surface(self, surface, delta_s=None):
		"""Broaden charge density isosurface."""
		if delta_s is None:
			delta_s = self.default_delta_s
		if abs(surface) < delta_s:
			return 15.0/(16.0*delta_s)*(1.0 - (surface/delta_s)**2)**2
		else:
			return 0.0

	def get_c(self, input_mesh, partial, fraction, delta_s):
		"""Return c mesh for broadened surface integration."""
		if partial:
			charge_density_mesh = abs(input_mesh)**2
		else:
			charge_density_mesh = input_mesh
		max_density = np.max(charge_density_mesh)
		isovalue = fraction * max_density

		log_mesh = np.empty(charge_density_mesh.shape, dtype=float)

		cd_zeros = charge_density_mesh == 0

		# Evaluate logarithm on unmasked entries
		log_mesh[~cd_zeros] = np.log(charge_density_mesh[charge_density_mesh != 0] / isovalue)
		log_mesh[cd_zeros] = np.inf

		# Apply broadening to surface
		broadened_mesh = np.where(abs(log_mesh) < delta_s, 15.0/(16.0*delta_s)*(1.0 - (log_mesh/delta_s)**2)**2, 0)

		# Remove overlapping layers of surface
		# Iterate over x and y
		for i, j in np.ndindex(broadened_mesh.shape[:2]):
			on_surface = False
			past_surface = False
			# Iterate over z, starting from top
			for k in reversed(range(broadened_mesh.shape[2])):
				if past_surface:
					# First surface has been traversed, so replace subsequent elements with zeros
					broadened_mesh[i, j, k] = 0
				elif broadened_mesh[i, j, k] != 0 and not on_surface:
					# Highest surface has been reached
					on_surface = True
				elif on_surface and broadened_mesh[i, j, k] == 0:
					# Was on surface, now just below
					on_surface = False
					past_surface = True

		# Get direction of gradient of isosurface
		gradient_surface = self.periodic_gradient(charge_density_mesh)

		# for i in range(3):
		# 	unit_vector_surface[i, gradient_magnitude_surface != 0] = gradient_surface[i, gradient_magnitude_surface != 0] / gradient_magnitude_surface[gradient_magnitude_surface != 0]
		vector_surface = np.zeros(gradient_surface.shape, dtype=float)

		for i in range(3):
			gradient_surface[~cd_zeros, i] = gradient_surface[~cd_zeros, i] / charge_density_mesh[~cd_zeros]
			gradient_surface[cd_zeros, i] = 0
			vector_surface[..., i] = np.multiply(broadened_mesh, gradient_surface[..., i])

		return vector_surface

	def get_A_mesh(self, c, wavefunction_mesh):
		grad_wavefunction = self.periodic_gradient(wavefunction_mesh)
		return self.mesh_dot_product(c, grad_wavefunction)

	@staticmethod
	def mesh_dot_product(vector_mesh_A, vector_mesh_B):
		"""Return dot/scalar product of two vector meshes"""
		# Check if resulting mesh will be complex
		if vector_mesh_A.dtype == complex or vector_mesh_B.dtype == complex:
			dtype = complex
		else:
			dtype = float

		scalar_mesh = np.zeros(vector_mesh_A.shape[:3], dtype=dtype)
		for i in range(3):
			scalar_mesh += vector_mesh_A[..., i]*vector_mesh_B[..., i]
		return scalar_mesh

	@staticmethod
	def grid_dot_product(vector_mesh_A, vector_mesh_B):
		"""Return dot/scalar product of two vector meshes"""
		# Check if resulting mesh will be complex
		if vector_mesh_A.dtype == complex or vector_mesh_B.dtype == complex:
			dtype = complex
		else:
			dtype = float

		scalar_mesh = np.zeros(vector_mesh_A.shape[:2], dtype=dtype)
		for i in range(3):
			scalar_mesh += vector_mesh_A[..., i]*vector_mesh_B[..., i]
		return scalar_mesh

	def get_B_mesh(self, c, wavefunction_mesh):
		B = np.zeros_like(c, dtype=complex)
		for i in range(3):
			B[..., i] = c[..., i] * wavefunction_mesh
		return B

	def greens_function_mesh(self, z_index, tip_work_func, tip_energy, debug=False):

		plane_length_half = - (- (max(self.real_mesh.shape[:2]) + 1) / 2)
		plane_length_full = plane_length_half * 2 - 1

		plane_shape_half = (plane_length_half,) * 2

		plane_UR = np.zeros(plane_shape_half, dtype=float)

		G_shape = (plane_length_full, plane_length_full, self.real_mesh.shape[2])
		G_mesh = np.zeros(G_shape, dtype=float)

		print self.real_mesh.shape, plane_length_half, plane_length_full, G_shape

		debug_str = "Calculating G(r - R): "
		if debug:
			sys.stdout.write(debug_str)
			sys.stdout.flush()
		points_done = 0
		bars_done = 0

		for k in range(G_shape[2]):
			for i in range(plane_length_half):
				for j in range(i + 1):
					if k >= z_index:
						G = 0
					else:
						distance = self.grid_spacing * np.sqrt(i**2 + j**2 + (z_index - k)**2)
						G = self.greens_function(distance, tip_work_func, tip_energy)

					for x, y in [i, j], [j, i]:
						plane_UR[x, y] = G

			plane_UL = np.flipud(plane_UR[1:, :])
			plane_LR = np.fliplr(plane_UR[:, 1:])
			plane_LL = np.flipud(plane_LR[1:, :])

			plane_U = np.concatenate((plane_UL, plane_UR), axis=0)
			plane_L = np.concatenate((plane_LL, plane_LR), axis=0)
			plane = np.concatenate((plane_L, plane_U), axis=1)

			G_mesh[..., k] = plane

			points_done += 1
			prog = float(points_done) / self.real_mesh.shape[2]
			if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
				percent = prog * 100
				sys.stdout.write('\r')
				sys.stdout.write(debug_str)
				sys.stdout.write(" [{:<{}}]".format(self.PROG_BAR_CHARACTER * bars_done, self.PROG_BAR_INTERVALS))
				sys.stdout.write(" {:3.0f}%".format(percent))
				sys.stdout.flush()
				bars_done += 1

		if debug:
			sys.stdout.write("\n")
			sys.stdout.flush()
		return G_mesh

	def propagated_psi_filename(self, K, E, T, fraction, delta_s):
		"""Return standardised filename for relevant propagated wavefunction file"""
		return (self.MESH_FOLDER+self.PROP_PSI_FNAME+self.name+"_"+str(self.grid_spacing)+"_"
				+str(K.x)+"_"+str(K.y)+"_"+str(K.z)+"_"+str(E)+"_"+str(T)+"_"+str(fraction)+"_"+str(delta_s)+self.EXT)

	def write_prop_psi(self, psi, K, E, T, fraction, delta_s):
		"""Write wavefunction function mesh to file"""
		filename = self.propagated_psi_filename(K, E, T, fraction, delta_s)
		psi_file = safe_open(filename, "w")
		for i, j in np.ndindex(psi.shape):
			if psi[i, j] != 0:
				p = psi[i, j]
				psi_file.write(str(i)+" "+str(j)+" "+str(p.real)+" "+str(p.imag)+"\n")
		psi_file.close()

	def read_prop_psi(self, K, E, T, fraction, delta_s=None, debug=False, debug_file=False):
		"""Read wavefunction mesh from file"""
		if delta_s is None:
			delta_s = self.default_delta_s

		filename = self.propagated_psi_filename(K, E, T, fraction, delta_s)
		psi_file = open(filename, "r")
		psi = np.zeros(self.real_mesh.shape[:2], dtype=complex)

		if debug_file:
			sys.stdout.write("Reading psi(R) from {}\n".format(filename))
			sys.stdout.flush()
		elif debug:
			if self.PRINT_RELATIVE_TO_EF:
				E_str = str(E - self.fermi_level) + " eV"
			else:
				E_str = str(E) + " eV"
			sys.stdout.write("Reading psi(R) at k = {!s}, E = {}\n".format(K, E_str))
			sys.stdout.flush()

		for line in psi_file:
			line_split = line.split()
			i = int(line_split[0])
			j = int(line_split[1])
			real = float(line_split[2])
			imag = float(line_split[3])
			psi[i, j] = complex(real, imag)
		return psi

	def calculate_current_scan_iso(self, z, V, T, tip_work_func, tip_energy, delta_s=None, fraction=0.025, recalculate=False, write=True, vectorised=True, partial_surface=False, debug=False):
		"""Calculate tunnelling current across plane.

		Args:
			z (float): z-value of plane
			V (float): Bias voltage
			T (float): Absolute temperature
			tip_work_func (float): Work function of tip
			tip_energy (float): Fermi-level of tip
			delta_s (float): Surface width
			fraction (float, opt.): Fraction of maximum charge density to use as isovalue for isosurface
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file
			debug (bool, opt.): Print extra information during runtime
		"""

		if delta_s is None:
			delta_s = self.default_delta_s

		if debug:
			sys.stdout.write("Calculating I(R)\n")
			sys.stdout.flush()

		if V > 0:
			min_E = self.fermi_level
			max_E = self.fermi_level + V
		else:
			min_E = self.fermi_level + V
			max_E = self.fermi_level

		total_k_weight = 0
		total_energies = 0
		for K in self.bands:
			total_k_weight += K.weight
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					total_energies += 1
		energies_done = 0

		z, k = self.get_nearest_mesh_value(z, indices=True, points=self.real_mesh.shape[2])
		current = np.zeros(self.real_mesh.shape[:2], dtype=float)
		psi = np.zeros_like(current, dtype=complex)
		elements = current.shape[0]*current.shape[1]

		if not partial_surface:
			ldos = self.get_ldos_grid(min_E, max_E, T, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)
			c = self.get_c(ldos, False, fraction, delta_s)
		G_conjugate = np.conjugate(self.greens_function_mesh(k, tip_work_func, tip_energy, debug=debug))

		if debug:
			sys.stdout.write("Calculating grad(G)\n")
			sys.stdout.flush()
		G_conjugate_gradient = self.periodic_gradient(G_conjugate)

		G_centre = G_conjugate.shape[0] / 2

		for K in self.bands:
			w = K.weight
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					fd = self.fermi_dirac(E, T)
					if V > 0:
						fd = 1 - fd

					if not recalculate and os.path.isfile(self.propagated_psi_filename(K, E, T, fraction, delta_s)):
						# Read data from file
						psi = self.read_prop_psi(K, E, T, fraction, delta_s=delta_s, debug=debug)
						read = True
					else:
						read = False
						points_done = 0
						bars_done = 0

						raw_psi = self.get_psi_grid(K, E, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)

						prog = float(energies_done) / total_energies * 100
						if self.PRINT_RELATIVE_TO_EF:
							E_str = str(E - self.fermi_level) + " eV"
						else:
							E_str = str(E) + " eV"
						debug_str = "Calculating psi(R) at k = {!s}, E = {}: {:5.1f}%".format(K, E_str, prog)
						if debug:
							sys.stdout.write(debug_str)
							sys.stdout.flush()
						if partial_surface:
							c = self.get_c(raw_psi, True, fraction, delta_s)

						A = self.get_A_mesh(c, raw_psi)
						B = self.get_B_mesh(c, raw_psi)

						for i, j in np.ndindex(current.shape):

							G_conjugate_rolled = np.roll(G_conjugate, (i - G_centre), 0)
							G_conjugate_rolled = np.roll(G_conjugate_rolled, (j - G_centre), 1)[
							                     :self.real_mesh.shape[0], :self.real_mesh.shape[1]]

							G_conjugate_gradient_rolled = np.roll(G_conjugate_gradient, (i - G_centre), 0)
							G_conjugate_gradient_rolled = np.roll(G_conjugate_gradient_rolled, (j - G_centre), 1)[
							                              :self.real_mesh.shape[0], :self.real_mesh.shape[1]]

							integrand = G_conjugate_rolled * A - self.mesh_dot_product(B, G_conjugate_gradient_rolled)
							psi[i, j] = np.sum(integrand) * self.grid_spacing ** 3

							points_done += 1
							prog = float(points_done) / elements

							if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
								percent = prog * 100
								sys.stdout.write('\r')
								sys.stdout.write(debug_str)
								sys.stdout.write(" [{:<{}}] {:3.0f}%".format(self.PROG_BAR_CHARACTER * bars_done,
								                                             self.PROG_BAR_INTERVALS, percent))
								sys.stdout.flush()
								bars_done += 1

						if write:
							self.write_prop_psi(psi, K, E, T, fraction, delta_s)

					current += fd * (w / total_k_weight) * abs(psi)**2

					energies_done += 1

					if debug and not read:
						sys.stdout.write("\n")
						sys.stdout.flush()
		return current

	def calculate_current_scan_flat(self, z, V, T, tip_work_func, tip_energy, wf_height, recalculate=False, write=True, vectorised=True, partial_surface=False, debug=False):
		"""Calculate tunnelling current across plane.

		Args:
			z (float): z-value of plane
			V (float): Bias voltage
			T (float): Absolute temperature
			tip_work_func (float): Work function of tip
			tip_energy (float): Fermi-level of tip
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file
			debug (bool, opt.): Print extra information during runtime
		"""

		if debug:
			sys.stdout.write("Calculating I(R) at V = {} eV\n".format(V))
			sys.stdout.flush()

		if V > 0:
			min_E = self.fermi_level
			max_E = self.fermi_level + V
		else:
			min_E = self.fermi_level + V
			max_E = self.fermi_level

		total_k_weight = 0
		total_energies = 0
		for K in self.bands:
			total_k_weight += K.weight
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					total_energies += 1
		if debug:
			sys.stdout.write("Total Energies = {}\n".format(total_energies))
			sys.stdout.flush()
		energies_done = 0

		z, k = self.get_nearest_mesh_value(z, indices=True, points=self.real_mesh.shape[2])
		wf_height, wf_index = self.get_nearest_mesh_value(wf_height, indices=True, points=self.real_mesh.shape[2])
		current = np.zeros(self.real_mesh.shape[:2], dtype=float)
		psi = np.zeros_like(current, dtype=complex)

		G_conjugate = np.conjugate(self.greens_function_mesh(k, tip_work_func, tip_energy, debug=debug))

		if debug:
			sys.stdout.write("Calculating grad(G)")
			sys.stdout.flush()
		G_conjugate_gradient = self.periodic_gradient(G_conjugate)

		G_conjugate = G_conjugate[..., wf_index]
		G_conjugate_gradient = G_conjugate_gradient[..., wf_index, :]
		G_centre = G_conjugate.shape[0] / 2

		for K in self.bands:
			w = K.weight
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					fd = self.fermi_dirac(E, T)
					if V > 0:
						fd = 1 - fd

					raw_psi = self.get_psi_grid(K, E, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)
					grad_raw_psi = self.periodic_gradient(raw_psi)[..., wf_index, :]
					raw_psi = raw_psi[..., wf_index]

					prog = float(energies_done) / total_energies * 100
					if self.PRINT_RELATIVE_TO_EF:
						E_str = str(E - self.fermi_level) + " eV"
					else:
						E_str = str(E) + " eV"
					debug_str = "Calculating psi(R) at k = {!s}, E = {}; {:5.1f}%\n".format(K, E_str, prog)
					if debug:
						sys.stdout.write(debug_str)
						sys.stdout.flush()

					for i, j in np.ndindex(current.shape):

						G_conjugate_rolled = np.roll(G_conjugate, (i - G_centre), 0)
						G_conjugate_rolled = np.roll(G_conjugate_rolled, (j - G_centre), 1)[
						                     :self.real_mesh.shape[0], :self.real_mesh.shape[1]]

						G_conjugate_gradient_rolled = np.roll(G_conjugate_gradient, (i - G_centre), 0)
						G_conjugate_gradient_rolled = np.roll(G_conjugate_gradient_rolled, (j - G_centre), 1)[
						                              :self.real_mesh.shape[0], :self.real_mesh.shape[1]]

						integrand = G_conjugate_rolled * raw_psi - self.grid_dot_product(grad_raw_psi, G_conjugate_gradient_rolled)
						psi[i, j] = np.sum(integrand) * self.grid_spacing ** 2

					current += fd * (w / total_k_weight) * abs(psi)**2

					energies_done += 1
		return current

	def current_filename(self, z, V, T):
		"""Return standardised filename for relevant current file"""
		return self.MESH_FOLDER+self.CURRENT_FNAME+self.name+"_"+str(self.grid_spacing)+"_"+str(z)+"_"+str(V)+"_"+str(T)+self.EXT

	def current_plane_filename(self, z, wf_height, V, T):
		"""Return standardised filename for relevant current file"""
		return self.MESH_FOLDER+self.CURRENT_FNAME+self.name+"_"+str(self.grid_spacing)+"_"+str(z)+"_"+str(wf_height)+"_"+str(V)+"_"+str(T)+self.EXT

	def write_current(self, current, z, V, T, debug=False):
		"""Write current to file.

		Args:
			min_E: Minimum energy
			max_E: Maximum energy
			T: Absolute temperature in K
			debug (bool, opt.): Print extra information during runtime
		"""
		filename = self.current_filename(z, V, T)

		current_file = safe_open(filename, "w")
		if debug:
			sys.stdout.write("Writing current grid to {}\n".format(filename))
		# Iterate over mesh points
		for i, j in np.ndindex(current.shape):
			# If LDOS is non-zero at mesh point, write data to file
			if current[i, j] != 0:
				current_file.write(str(i)+" "+str(j)+" "+str(current[i, j])+"\n")
		current_file.close()

	def write_current_plane(self, current, z, wf_height, V, T, debug=False):
		"""Write current to file.

		Args:
			T: Absolute temperature in K
			debug (bool, opt.): Print extra information during runtime
		"""
		filename = self.current_plane_filename(z, wf_height, V, T)

		current_file = safe_open(filename, "w")
		if debug:
			sys.stdout.write("Writing current grid to {}\n".format(filename))
		# Iterate over mesh points
		for i, j in np.ndindex(current.shape):
			# If LDOS is non-zero at mesh point, write data to file
			if current[i, j] != 0:
				current_file.write(str(i)+" "+str(j)+" "+str(current[i, j])+"\n")
		current_file.close()

	def read_current(self, z, V, T, debug=False):
		"""Read current grid from file"""
		filename = self.current_filename(z, V, T)
		current_file = open(filename, 'r')
		current = np.zeros(self.real_mesh.shape[:2], dtype=float)

		if debug:
			sys.stdout.write("Reading I(R) from {}\n".format(filename))
			sys.stdout.flush()

		for line in current_file:
			line_split = line.split()
			# Get mesh indices
			i = int(line_split[0])
			j = int(line_split[1])

			# Get current value
			value = float(line_split[2])
			current[i, j] = value

		if debug:
			sys.stdout.write("I(R) successfully read\n")
			sys.stdout.flush()
		return current

	def read_current_plane(self, z, wf_index, V, T, debug=False):
		"""Read current grid from file"""
		filename = self.current_plane_filename(z, wf_index, V, T)
		current_file = open(filename, 'r')
		current = np.zeros(self.real_mesh.shape[:2], dtype=float)

		if debug:
			sys.stdout.write("Reading I(R) from {}\n".format(filename))
			sys.stdout.flush()

		for line in current_file:
			line_split = line.split()
			# Get mesh indices
			i = int(line_split[0])
			j = int(line_split[1])

			# Get current value
			value = float(line_split[2])
			current[i, j] = value

		if debug:
			sys.stdout.write("I(R) successfully read\n")
			sys.stdout.flush()
		return current

	def get_current_scan_iso(self, z, V, T, tip_work_func, tip_energy, delta_s=None, fraction=0.025, recalculate=False, write=True, vectorised=True, partial_surface=False, debug=False):
		if delta_s is None:
			delta_s = self.default_delta_s
		if not recalculate and os.path.isfile(self.current_filename(z, V, T)):
			# Read data from file
			current = self.read_current(z, V, T, debug=debug)
		else:
			current = self.calculate_current_scan_iso(z, V, T, tip_work_func, tip_energy, delta_s=delta_s, fraction=fraction, recalculate=recalculate, write=write, vectorised=vectorised, partial_surface=partial_surface, debug=debug)
			if write:
				self.write_current(current, z, V, T)
		return current

	def get_current_scan_plane(self, z, wf_height, V, T, tip_work_func, tip_energy, recalculate=False, write=True, vectorised=True, partial_surface=False, debug=False):
		if not recalculate and os.path.isfile(self.current_plane_filename(z, wf_height, V, T)):
			# Read data from file
			current = self.read_current_plane(z, wf_height, V, T, debug=debug)
		else:
			current = self.calculate_current_scan_flat(z, V, T, tip_work_func, tip_energy, wf_height, recalculate=recalculate, write=write, vectorised=vectorised, partial_surface=partial_surface, debug=debug)
			if write:
				self.write_current_plane(current, z, wf_height, V, T)
		return current

	def get_spectrum(self, x, y, z, min_V, max_V, T, dE=0.005, debug=False):

		min_E = min_V + self.fermi_level
		max_E = max_V + self.fermi_level

		x, i = self.get_nearest_mesh_value(x, indices=True, points=self.real_mesh.shape[0])
		y, j = self.get_nearest_mesh_value(y, indices=True, points=self.real_mesh.shape[1])
		z, k = self.get_nearest_mesh_value(z, indices=True, points=self.real_mesh.shape[2])

		sigma = self.BOLTZMANN * T

		Es = []
		psis = []
		weights = []
		for K in self.bands:
			for E in self.bands[K]:
				if min_E - 3*sigma < E < max_E + 3*sigma:
					Es.append(E)
					weights.append(K.weight)
					psi = self.get_psi_grid(K, E, debug=debug)[i, j, k]
					psis.append(abs(psi)**2)

		Es = np.array(Es)
		psis = np.array(psis)
		weights = np.array(weights)

		E_range = np.arange(min_E, max_E, dE)
		LDOS = np.zeros_like(E_range)

		for u in range(len(E_range)):
			LDOS[u] = np.sum(weights * psis * np.exp(- (((E_range[u] - Es) / sigma)**2) / 2))
		V_range = E_range - self.fermi_level
		return V_range, LDOS

	def get_line_cut(self, axis, value, z, min_V, max_V, T, dE=0.0005, debug=False):

		if axis == 'x':
			i = np.arange(self.real_mesh.shape[0])
			y, j = self.get_nearest_mesh_value(value, indices=True, points=self.real_mesh.shape[1])
		elif axis == 'y':
			j = np.arange(self.real_mesh.shape[1])
			x, i = self.get_nearest_mesh_value(value, indices=True, points=self.real_mesh.shape[0])
		else:
			raise ValueError('Axis must be x or y')
		z, k = self.get_nearest_mesh_value(z, indices=True, points=self.real_mesh.shape[2])

		min_E = min_V + self.fermi_level
		max_E = max_V + self.fermi_level

		sigma = self.BOLTZMANN * T

		Es = []
		psis = []
		weights = []
		for K in self.bands:
			for E in self.bands[K]:
				if min_E - 3 * sigma < E < max_E + 3 * sigma:
					Es.append(E)
					weights.append(K.weight)
					psi = self.get_psi_grid(K, E, debug=debug)[i, j, k]
					psis.append(abs(psi) ** 2)

		Es = np.array(Es)
		psis = np.array(psis)
		print psis.shape
		weights = np.array(weights)

		E_range = np.arange(min_E, max_E, dE)
		LDOS = np.zeros((E_range.shape[0], psis.shape[1]))

		for u in range(len(E_range)):
			for v in range(LDOS.shape[1]):
				LDOS[u, v] = np.sum(weights * psis[:, v] * np.exp(- (((E_range[u] - Es) / sigma)**2) / 2))

		V_range = E_range - self.fermi_level

		return V_range, LDOS
