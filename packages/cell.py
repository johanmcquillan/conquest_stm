
import os
import sys
import numpy as np
import math
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
		x_mesh (3D mgrid): Mesh of x values
		y_mesh (3D mgrid): Mesh of y values
		z_mesh (3D mgrid): Mesh of z values
		atoms (int : Atom): Atom objects of simulation indexed by atom number
		bands (Vector : [float]): Energies of bands indexed by k-point vector
		support_grid (array(SmartDict)): Mesh of support function values, indexed by [x, y, z][atom_key][l][zeta][m]
	"""

	BOLTZMANN = 8.6173303E-5  # Boltzmann's Constant in eV/K
	ELECTRON_MASS = 9.10938E-31  # Electron Mass in kg
	H_BAR = 4.135667662E-15  # Reduced Planck's Constant in eV.s

	MESH_FOLDER = "meshes/"
	SUPPORT_FNAME = "supp_"
	LDOS_FNAME = "ldos_"
	PSI_FNAME = "psi_"
	CURRENT_FNAME = "current_"
	EXT = ".dat"

	PRINT_RELATIVE_TO_EF = True
	PROG_BAR_INTERVALS = 20
	PROG_BAR_CHARACTER = ">"

	def __init__(self, name, fermi_level, x_length, y_length, z_length, grid_spacing=0.5):
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

		vector_x = int(x_length / grid_spacing) * grid_spacing
		vector_y = int(y_length / grid_spacing) * grid_spacing
		vector_z = int(z_length / grid_spacing) * grid_spacing
		self.vector = Vector(vector_x, vector_y, vector_z)

		# Form Cartesian meshes
		self.real_mesh = np.transpose(np.mgrid[0: x_length: grid_spacing, 0: y_length: grid_spacing, 0: z_length: grid_spacing], (1, 2, 3, 0))
		self.vector_mesh = self.get_vector_mesh()
		self.mesh_points = self.real_mesh.shape[0] * self.real_mesh.shape[1] * self.real_mesh.shape[2]

		# Initialise atoms and bands
		self.atoms = {}
		self.bands = {}
		self.support_grid = None

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
		"""Calculate Fermi-Dirac distribution value.

		Args:
			energy (float): Energy in eV
			temperature (float): Absolute temperature in K

		Returns:
			float: Occupation value
		"""
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

	def calculate_support_grid(self, vectorised=True, debug=False):
		"""Evaluate support function for each PAO on mesh.

		Args:
			debug (bool, opt.): If true, print extra information during runtime

		Returns:
			array(SmartDict): Mesh of support function values, indexed by [x, y, z][atom_key][l][zeta][m]
		"""
		if debug:
			print "Calculating support grid"
		# Initialise support grid
		support_grid = np.empty_like(self.real_mesh[..., 0], dtype=SmartDict)

		# Iterate over all atoms
		for atom_key in self.atoms:
			atom = self.atoms[atom_key]

			debug_str = "    Support grid for atom {:3} - {!s}: ".format(atom_key, atom.ion_name)
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

			# Iterate over mesh points within cutoff
			i_end = i_start + lm_shape[0] + 2
			j_end = j_start + lm_shape[1] + 2
			k_end = k_start + lm_shape[2] + 2

			rolled_mesh = np.roll(support_grid, -i_start, 0)
			rolled_mesh = np.roll(rolled_mesh, -j_start, 1)
			rolled_mesh = np.roll(rolled_mesh, -k_start, 2)
			partial_mesh = rolled_mesh[0:lm_shape[0], 0:lm_shape[1], 0:lm_shape[2]]
			if not False:
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
			else:
				rolled_mesh = np.roll(support_grid, -i_start, 0)
				rolled_mesh = np.roll(rolled_mesh, -j_start, 1)
				rolled_mesh = np.roll(rolled_mesh, -k_start, 2)

				partial_mesh = rolled_mesh[0:lm_shape[0], 0:lm_shape[1], 0:lm_shape[2]]
				vector_mesh = Vector.array(local_mesh) - atom.atom_pos
				print vector_mesh.shape, partial_mesh.shape, self.real_mesh.shape
				pg = atom.calculate_support_point_vec(vector_mesh, partial_mesh, atom_key)

				rolled_mesh[0:lm_shape[0], 0:lm_shape[1], 0:lm_shape[2]] = partial_mesh
				rolled_mesh = np.roll(rolled_mesh, i_start, 0)
				rolled_mesh = np.roll(rolled_mesh, j_start, 1)
				support_grid = np.roll(rolled_mesh, k_start, 2)

				print "Calculated "+str(atom_key)

		return support_grid

	def support_filename(self):
		"""Return standardised filename for relevant support function file"""
		return self.MESH_FOLDER+self.SUPPORT_FNAME+self.name+"_"+str(self.grid_spacing)+self.EXT

	def write_support_grid(self, debug=False):
		"""Write support function mesh to file"""
		filename = self.support_filename()
		support_file = safe_open(filename, 'w')

		if debug:
			sys.stdout.write("Writing support grid to "+filename+": ")
			sys.stdout.flush()
		points_done = 0
		bars_done = 0

		# Iterate over mesh points
		for ijk in np.ndindex(self.real_mesh.shape[:3]):
			i, j, k = ijk
			# If support function values exist at mesh point
			if self.support_grid[ijk]:
				support_file.write(str(i) + " " + str(j) + " " + str(k) + "\n")
				# Iterate over atoms
				for atom_key in self.support_grid[ijk]:
					# Write atom index
					support_file.write(str(atom_key) + "\n")
					# Iterate over orbitals
					for l in self.support_grid[ijk][atom_key]:
						for zeta in self.support_grid[ijk][atom_key][l]:
							for m in self.support_grid[ijk][atom_key][l][zeta]:
								# Write orbital data
								line = (str(l) + " " + str(zeta) + " " + str(m) + " "
										+ str(self.support_grid[ijk][atom_key][l][zeta][m]))
								support_file.write(line + "\n")

			points_done += 1
			if debug and float(points_done) / self.mesh_points * self.PROG_BAR_INTERVALS > bars_done:
				sys.stdout.write(self.PROG_BAR_CHARACTER)
				sys.stdout.flush()
				bars_done += 1
		if debug:
			print
		support_file.close()

	def read_support_grid(self, debug=False):
		"""Read support function mesh from file"""
		filename = self.support_filename()
		support_file = open(filename, 'r')
		support_grid = np.empty_like(self.real_mesh[..., 0], dtype=SmartDict)

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

	def get_support_grid(self, recalculate=False, write=True, vectorised=True, debug=False):
		"""Get support function mesh.

		Args:
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file
			debug (bool, opt.): Print extra information during runtime

		Returns:
			array(SmartDict): Mesh of support function values, indexed by [x, y, z][atom_key][l][zeta][m]
		"""
		if self.support_grid is None:
			if not recalculate and os.path.isfile(self.support_filename()):
				# Read support grid from file
				self.support_grid = self.read_support_grid(debug=debug)
			else:
				# Recalculate support grid
				self.support_grid = self.calculate_support_grid(vectorised=vectorised, debug=debug)
				# Write to file
				if write:
					self.write_support_grid(debug=debug)
		return self.support_grid

	def calculate_psi_grid_vec(self, support_dict, atom_key, l, zeta, m, coefficient):
		# Evaluate wavefunction contribution
		if support_dict and atom_key in support_dict:
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
		support_grid = self.get_support_grid(recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)

		if self.PRINT_RELATIVE_TO_EF:
			E_str = str(E - self.fermi_level) + " eV"
		else:
			E_str = str(E) + " eV"
		debug_str = "Calculating psi(r) at k = "+str(K)+", E = "+E_str+": "

		if debug:
			sys.stdout.write(debug_str)
			sys.stdout.flush()

		for atom_key in self.atoms:
			atom = self.atoms[atom_key]
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
			print
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
			print "Reading psi(r) from "+filename
		elif debug:
			if self.PRINT_RELATIVE_TO_EF:
				E_str = str(E - self.fermi_level) + " eV"
			else:
				E_str = str(E) + " eV"
			print "Reading psi(r) at k = "+str(K)+", E = "+E_str

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
			debug (bool, opt.): Print extra information during runtime

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

	def calculate_ldos_grid(self, min_E, max_E, T, recalculate=False, write=True, vectorised=False, debug=False):
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
			print "Calculating local density of states grid"
		ldos_grid = np.zeros_like(self.real_mesh[..., 0], dtype=float)

		total_k_weight = 0
		for K in self.bands:
			total_k_weight += K.weight

		for K in self.bands:
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					psi_grid = self.get_psi_grid(K, E, recalculate=recalculate, write=write, debug=debug)
					fd = self.fermi_dirac(E, T)
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
		for ijk in np.ndindex(self.real_mesh.shape[:3]):
			# If LDoS is non-zero at mesh point, write data to file
			i, j, k = ijk
			if ldos_grid[ijk]:
				ldos_file.write(str(i)+" "+str(j)+" "+str(k)+" "+str(ldos_grid[i, j, k])+"\n")
		ldos_file.close()

	def read_ldos_grid(self, min_E, max_E, T, debug=False):
		"""Read LDoS mesh from file"""
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
		"""Get LDoS mesh.

		Args:
			min_E: Minimum energy
			max_E: Maximum energy
			T: Absolute temperature in K
			recalculate (bool, opt.): Force recalculation of meshes, even if already stored
			write (bool, opt.): Write calculated grids to file
			debug (bool, opt.): Print extra information during runtime

		Returns:
			array(float): Mesh of LDOS values
		"""
		# Read ldos grid from file if not stored by cell
		if not recalculate and os.path.isfile(self.ldos_filename(min_E, max_E, T)):
			ldos_grid = self.read_ldos_grid(min_E, max_E, T, debug=debug)
		else:
			# Calculate LDoS on mesh
			ldos_grid = self.calculate_ldos_grid(min_E, max_E, T, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)
			if write:
				self.write_ldos_grid(ldos_grid, min_E, max_E, T, debug=debug)
		return ldos_grid

	def get_vector_mesh(self):
		vector_mesh = np.empty_like(self.real_mesh[..., 0], dtype=Vector)
		for ijk in np.ndindex(self.real_mesh.shape[:3]):
			x, y, z = self.real_mesh[ijk]
			r = Vector(x, y, z)
			vector_mesh[ijk] = r
		return vector_mesh

	def kappa_squared(self, tip_work_func, tip_fermi_level):
		return 2*self.ELECTRON_MASS/(self.H_BAR**2)*(tip_work_func - tip_fermi_level)

	def greens_function(self, distance, tip_work_func, tip_energy):
		if distance == 0:
			return 0
		else:
			kappa2 = self.kappa_squared(tip_work_func, tip_energy)
			return np.exp(- kappa2 * distance) / (4*np.pi*distance)

	def filtered_greens_function(self, position, tip_position, tip_work_func, tip_energy, alpha):
		distance = abs(position - tip_position)
		kappa2 = self.kappa_squared(tip_work_func, tip_energy)
		kappa = np.sqrt(kappa2)
		u = np.sqrt(2*np.pi)/(4*distance)*np.exp(alpha**2*kappa2)
		v = np.exp(np.sqrt(kappa)*distance)*math.erf(alpha*kappa + distance/(2*alpha))
		w = np.exp(- kappa * distance)*math.erf(alpha*kappa - distance/(2*alpha))
		x = 2*np.sinh(kappa * distance)
		return u*(v - w - x)

	def broadened_surface(self, surface, delta_s):
		if abs(surface) < delta_s:
			return 15.0/(16.0*delta_s)*(1.0 - (surface/delta_s)**2)**2
		else:
			return 0.0

	def get_c(self, wavefunction_mesh, fraction, delta_s):
		charge_density_mesh = abs(wavefunction_mesh)**2
		max_density = np.max(charge_density_mesh)
		isovalue = fraction * max_density

		masked_mesh = np.ma.masked_array(charge_density_mesh, mask=np.where(charge_density_mesh == 0, 1, 0))
		log_mesh = np.empty_like(charge_density_mesh)
		log_mesh[~masked_mesh.mask] = np.log(masked_mesh[~masked_mesh.mask] / isovalue)

		broadened_mesh = np.where(abs(masked_mesh) < delta_s, 15.0/(16.0*delta_s)*(1.0 - (masked_mesh/delta_s)**2)**2, 0)

		for ij in np.ndindex(broadened_mesh.shape[:2]):
			i, j = ij
			on_surface = False
			past_surface = False
			for k in reversed(range(broadened_mesh.shape[2])):
				if past_surface:
					broadened_mesh[i, j, k] = 0
				elif broadened_mesh[i, j, k] != 0 and not on_surface:
					on_surface = True
				elif on_surface and broadened_mesh[i, j, k] == 0:
					on_surface = False
					past_surface = True

		unit_vector_surface = np.array(np.gradient(charge_density_mesh, self.grid_spacing))
		vector_surface = np.transpose(np.multiply(broadened_mesh, unit_vector_surface), (1, 2, 3, 0))

		return vector_surface

	def get_A_mesh(self, c, wavefunction_mesh):
		grad_wavefunction = np.transpose(np.array(np.gradient(wavefunction_mesh, self.grid_spacing)), (1, 2, 3, 0))
		return self.mesh_dot_product(c, grad_wavefunction)

	@staticmethod
	def mesh_dot_product(vector_mesh_A, vector_mesh_B):
		if vector_mesh_A.dtype == complex or vector_mesh_B.dtype == complex:
			dtype = complex
		else:
			dtype = float

		scalar_mesh = np.zeros(vector_mesh_A.shape[:3], dtype=dtype)
		for i in range(3):
			scalar_mesh += vector_mesh_A[..., i]*vector_mesh_B[..., i]
		return scalar_mesh

	def get_B_mesh(self, c, wavefunction_mesh):
		B = np.zeros_like(c, dtype=complex)
		for i in range(3):
			B[..., i] = c[..., i] * wavefunction_mesh
		return B

	def differences_mesh(self, z_index, debug=False):

		diff_shape = self.real_mesh.shape[:2] + self.real_mesh.shape[:3]
		differences_mesh = np.zeros(diff_shape, dtype=float)

		plane_length_h = max(self.real_mesh.shape[:2])
		plane_length_f = plane_length_h*2
		plane_shape = (plane_length_f,)*2

		plane = np.zeros(plane_shape, dtype=float)

		if debug:
			sys.stdout.write("Calculating difference mesh: ")
			sys.stdout.flush()
		points_done = 0
		bars_done = 0

		for k in range(self.real_mesh.shape[2]):
			for i in range(plane_length_h):
				for j in range(i + 1):
					if k < z_index:
						distance = self.grid_spacing * np.sqrt(i**2 + j**2 + k**2)
					else:
						distance = 0
					plane[plane_length_h + i, plane_length_h + j] = distance
					plane[plane_length_h - i, plane_length_h - j] = distance
					plane[plane_length_h + i, plane_length_h - j] = distance
					plane[plane_length_h - i, plane_length_h + j] = distance
					plane[plane_length_h + j, plane_length_h + i] = distance
					plane[plane_length_h - j, plane_length_h - i] = distance
					plane[plane_length_h + j, plane_length_h - i] = distance
					plane[plane_length_h - j, plane_length_h + i] = distance

				for ij in np.ndindex(differences_mesh.shape[:2]):
					i, j = ij
					i_start = plane_length_h - i
					j_start = plane_length_h - j
					i_end = plane_length_f - i
					j_end = plane_length_f - j
					differences_mesh[i, j, ..., k] = plane[i_start:i_end, j_start:j_end]

			points_done += 1
			if debug and float(points_done) / self.real_mesh.shape[2] * self.PROG_BAR_INTERVALS > bars_done:
				sys.stdout.write(self.PROG_BAR_CHARACTER)
				sys.stdout.flush()
				bars_done += 1

		if debug:
			print
		return differences_mesh

	def greens_function_mesh_partial(self, z_index, tip_work_func, tip_energy, c, debug=False):

		G_shape = self.real_mesh.shape[:2] + self.real_mesh.shape[:3]
		G_mesh = np.zeros(G_shape, dtype=float)

		elements = G_shape[0]*G_shape[1]

		debug_str = "Calculating G(r - R): "
		if debug:
			sys.stdout.write(debug_str)
			sys.stdout.flush()
		points_done = 0
		bars_done = 0

		for ab in np.ndindex(G_shape[:2]):
			for ijk in np.ndindex(G_shape[2:]):
				if np.any(c[ijk]):
					distance = abs(self.vector_mesh[ab + (z_index,)] - self.vector_mesh[ijk])
					G_mesh[ab + ijk] = self.greens_function(distance, tip_work_func, tip_energy)

			points_done += 1
			prog = float(points_done) / elements
			if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
				percent = prog * 100
				sys.stdout.write('\r')
				sys.stdout.write(debug_str)
				sys.stdout.write(" [{:<{}}]".format(self.PROG_BAR_CHARACTER * bars_done, self.PROG_BAR_INTERVALS))
				sys.stdout.write(" {:3.0f}%".format(percent))
				sys.stdout.flush()
				bars_done += 1

		return G_mesh

	def greens_function_mesh_full(self, z_index, tip_work_func, tip_energy, debug=False):

		G_shape = self.real_mesh.shape[:2] + self.real_mesh.shape[:3]
		G_mesh = np.zeros(G_shape, dtype=float)

		plane_length_h = max(self.real_mesh.shape[:2])
		plane_length_f = plane_length_h*2
		plane_shape = (plane_length_f,)*2

		plane = np.zeros(plane_shape, dtype=float)

		debug_str = "Calculating G(r - R): "
		if debug:
			sys.stdout.write(debug_str)
			sys.stdout.flush()
		points_done = 0
		bars_done = 0

		for k in range(self.real_mesh.shape[2]):
			for i in range(plane_length_h):
				for j in range(i + 1):
					if k < z_index:
						distance = self.grid_spacing * np.sqrt(i**2 + j**2 + (z_index - k)**2)
						G = self.greens_function(distance, tip_work_func, tip_energy)
					else:
						G = 0

					for a in [-1, 1]:
						for b in [-1, 1]:
							plane[plane_length_h + a*i, plane_length_h + b*j] = G

					# plane[plane_length_h + i, plane_length_h + j] = G
					# plane[plane_length_h - i, plane_length_h - j] = G
					# plane[plane_length_h + i, plane_length_h - j] = G
					# plane[plane_length_h - i, plane_length_h + j] = G
					# plane[plane_length_h + j, plane_length_h + i] = G
					# plane[plane_length_h - j, plane_length_h - i] = G
					# plane[plane_length_h + j, plane_length_h - i] = G
					# plane[plane_length_h - j, plane_length_h + i] = G

				for ij in np.ndindex(G_mesh.shape[:2]):
					i, j = ij
					i_start = plane_length_h - i
					j_start = plane_length_h - j
					i_end = plane_length_f - i
					j_end = plane_length_f - j
					G_mesh[i, j, ..., k] = plane[i_start:i_end, j_start:j_end]

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

	def greens_function_mesh(self, z_index, tip_work_func, tip_energy, debug=False):

		plane_length_full = max(self.real_mesh.shape[:2])
		plane_length_half = - (- (max(self.real_mesh.shape[:2]) + 1) / 2)

		plane_shape_half = (plane_length_half,) * 2

		plane_UR = np.zeros(plane_shape_half, dtype=float)

		G_shape = (plane_length_full + 1, plane_length_full + 1, self.real_mesh.shape[2])
		G_mesh = np.zeros(G_shape, dtype=float)

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

	def calculate_current_scan(self, z, V, T, tip_work_func, tip_energy, delta_s, fraction=0.025, recalculate=False, write=True, vectorised=True, debug=False):
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
		elements = current.shape[0]*current.shape[1]

		# if debug:
		# 	sys.stdout.write("Calculating difference mesh: ")
		# 	sys.stdout.flush()
		# points_done = 0
		# bars_done = 0

		# print c.shape
		#
		# difference_dict = {}
		#
		# diff_shape = scan_vectors.shape + self.vector_mesh.shape
		# difference_array = np.zeros(diff_shape, dtype=float)
		# for indices_1 in np.ndindex(self.vector_mesh.shape):
		#
		# 	if indices_1[2] > self.vector_mesh.shape[2]/2:
		# 		r = self.vector_mesh[indices_1]
		#
		# 		for indices_2 in np.ndindex(scan_vectors.shape):
		# 			R = scan_vectors[indices_2]
		# 			difference = r - R
		# 			if difference in difference_dict:
		# 				distance = difference_dict[difference]
		# 			else:
		# 				distance = abs(difference)
		# 			difference_array[indices_2 + indices_1] = distance
		#
		# 		points_done += 1
		# 		if debug and float(points_done) / self.mesh_points * self.PROG_BAR_INTERVALS > bars_done:
		# 			sys.stdout.write(self.PROG_BAR_CHARACTER)
		# 			sys.stdout.flush()
		# 			bars_done += 1
		# if debug:
		# 	print

		# difference_mesh = np.zeros_like(self.vector_mesh, dtype=float)
		#
		# for i in range(self.vector_mesh.shape[2]):
		# 	vector_slice = self.vector_mesh[:, :, i]
		# 	difference_slice = abs(vector_slice - scan_vectors)
		# 	difference_mesh[:, :, i] = difference_slice

		# difference_array = self.differences_mesh(k, debug=debug)

		ldos = self.get_ldos_grid(min_E, max_E, T, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)
		c = self.get_c(ldos, fraction, delta_s)

		G_conjugate = np.conjugate(self.greens_function_mesh(k, tip_work_func, tip_energy, debug=debug))

		if debug:
			print "Calculating grad(G)"
		G_conjugate_gradient = np.transpose(np.array(np.gradient(G_conjugate, self.grid_spacing)), (1, 2, 3, 0))

		G_centre = G_conjugate.shape[0] / 2

		for K in self.bands:
			w = K.weight
			for E in self.bands[K]:
				if min_E <= E <= max_E:
					fd = self.fermi_dirac(E + V, T)
					wavefunction = self.get_psi_grid(K, E, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)

					prog = float(energies_done) / total_energies * 100
					if self.PRINT_RELATIVE_TO_EF:
						E_str = str(E - self.fermi_level) + " eV"
					else:
						E_str = str(E) + " eV"
					debug_str = "Calculating psi(R) at k = {!s}, E = {}: {:5.1f}%".format(K, E_str, prog)
					if debug:
						sys.stdout.write(debug_str)
						sys.stdout.flush()

					points_done = 0
					bars_done = 0

					A = self.get_A_mesh(c, wavefunction)
					B = self.get_B_mesh(c, wavefunction)

					for i, j in np.ndindex(current.shape):

						G_conjugate_rolled = np.roll(G_conjugate, (i - G_centre), 0)
						G_conjugate_rolled = np.roll(G_conjugate_rolled, (j - G_centre), 1)[:self.real_mesh.shape[0], :self.real_mesh.shape[1]]

						G_conjugate_gradient_rolled = np.roll(G_conjugate_gradient, (i - G_centre), 0)
						G_conjugate_gradient_rolled = np.roll(G_conjugate_gradient_rolled, (j - G_centre), 1)[:self.real_mesh.shape[0], :self.real_mesh.shape[1]]

						integrand = G_conjugate_rolled * A - self.mesh_dot_product(B, G_conjugate_gradient_rolled)
						psi = np.sum(integrand) * self.grid_spacing**3

						current[i, j] += fd * (w / total_k_weight) * abs(psi)**2
						points_done += 1
						prog = float(points_done) / elements

						if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
							percent = prog * 100
							sys.stdout.write('\r')
							sys.stdout.write(debug_str)
							sys.stdout.write(" [{:<{}}] {:3.0f}%".format(self.PROG_BAR_CHARACTER * bars_done, self.PROG_BAR_INTERVALS, percent))
							sys.stdout.flush()
							bars_done += 1

					energies_done += 1

					if debug:
						sys.stdout.write("\n")
						sys.stdout.flush()
		return current

	def current_filename(self, min_E, max_E, T):
		"""Return standardised filename for relevant LDOS file"""
		return self.MESH_FOLDER+self.CURRENT_FNAME+self.name+"_"+str(self.grid_spacing)+"_"+str(min_E)+"_"+str(max_E)+"_"+str(T)+self.EXT

	def write_current(self, current, min_E, max_E, T, debug=False):
		"""Write current to file.

		Args:
			min_E: Minimum energy
			max_E: Maximum energy
			T: Absolute temperature in K
			debug (bool, opt.): Print extra information during runtime
		"""
		filename = self.current_filename(min_E, max_E, T)

		current_file = safe_open(filename, "w")
		if debug:
			print "Writing current grid to "+filename
		# Iterate over mesh points
		for ij in np.ndindex(current.shape):
			# If LDoS is non-zero at mesh point, write data to file
			i, j = ij
			if current[ij] != 0:
				current_file.write(str(i)+" "+str(j)+" "+str(current[ij])+"\n")
		current_file.close()

	def read_current(self, min_E, max_E, T, debug=False):
		"""Read current mesh from file"""
		filename = self.current_filename(min_E, max_E, T)
		current_file = open(filename, 'r')
		current = np.zeros(self.real_mesh.shape[:2], dtype=float)

		if debug:
			print "Reading I(R) from "+filename

		for line in current_file:
			line_split = line.split()
			# Get mesh indices
			i = int(line_split[0])
			j = int(line_split[1])

			# Get current value
			value = float(line_split[2])
			current[i, j] = value

		if debug:
			print "current grid successfully read"
		return current

	def get_current_scan(self, z, V, T, tip_work_func, tip_energy, delta_s, fraction=0.025, recalculate=False, write=True, vectorised=True, debug=False):

		if V > 0:
			min_E = self.fermi_level
			max_E = self.fermi_level + V
		else:
			min_E = self.fermi_level + V
			max_E = self.fermi_level

		if not recalculate and os.path.isfile(self.current_filename(min_E, max_E, T)):
			# Read data from file
			current = self.read_current(min_E, max_E, T, debug=debug)
		else:
			current = self.calculate_current_scan(z, V, T, tip_work_func, tip_energy, delta_s, fraction=fraction, recalculate=recalculate, write=write, vectorised=vectorised, debug=debug)
			if write:
				self.write_current(current, min_E, max_E, T)
		return current

