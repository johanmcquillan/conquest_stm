
import os
import numpy as np
from copy import deepcopy

from sph import sph
from smartDict import SmartDict
from vector import Vector

BOLTZMANN = 8.6173303E-5  # Boltzmann's Constant in eV/K
SUPPORT_FNAME = "supp.aux"
LDOS_FNAME = "ldos.aux"


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
			xPoints (float): Number of points on x mesh
			yPoints (float): Number of points on y mesh
			zPoints (float): Number of points on z mesh
			xMesh (3D mgrid): Mesh of x values
			yMesh (3D mgrid): Mesh of y values
			zMesh (3D mgrid): Mesh of z values
			atoms (int : Atom): Atom objects of simulation indexed by atom number
			basisPoints (SmartDict): Basis function values, indexed by [atomKey][x][y][z][l][zeta][m]
			bands (Vector : [float]): Energies of bands indexed by k-point vector
	"""

	def __init__(self, name, fermiLevel, xLength, yLength, zLength, grid_spacing=0.5, debug=False):
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

		# Calculator number of points
		self.xPoints = int(self.vector.x / grid_spacing)
		self.yPoints = int(self.vector.y / grid_spacing)
		self.zPoints = int(self.vector.z / grid_spacing)

		# Form Cartesian meshes
		self.xMesh, self.yMesh, self.zMesh = np.mgrid[0: xLength: grid_spacing, 0: yLength: grid_spacing, 0: zLength: grid_spacing]

		# Initialise atoms and bands
		self.atoms = {}
		self.basisPoints = SmartDict()
		self.bands = {}
		self.support_grid = None
		self.ldos_grid = None

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

	def constrain_relative_vector(self, vector):
		"""Return a vector that is constrained within simulation cell"""
		new_vector = deepcopy(vector)

		# Check if vector components are greater than half of cell sides
		# If greater, add or subtract cell length

		if vector.x > self.vector.x/2:
			new_vector -= self.vector.project_x()
		elif vector.x <= -self.vector.x/2:
			new_vector += self.vector.project_x()

		if vector.y > self.vector.y/2:
			new_vector -= self.vector.project_y()
		elif vector.y <= -self.vector.y/2:
			new_vector += self.vector.project_y()

		if vector.z > self.vector.z/2:
			new_vector -= self.vector.project_z()
		elif vector.z <= -self.vector.z/2:
			new_vector += self.vector.project_z()

		# If vector is unchanged return
		# If vector has changed, constrain
		if new_vector == vector:
			return new_vector
		else:
			return self.constrain_relative_vector(new_vector)


	def constrain_vector_to_cell(self, vector):
		"""Return a vector that is constrained within simulation cell"""
		new_vector = deepcopy(vector)

		# Check if vector components are greater than half of cell sides
		# If greater, add or subtract cell length

		if vector.x >= self.vector.x:
			new_vector -= self.vector.project_x()
		elif vector.x < 0:
			new_vector += self.vector.project_x()

		if vector.y >= self.vector.y:
			new_vector -= self.vector.project_y()
		elif vector.y < 0:
			new_vector += self.vector.project_y()

		if vector.z >= self.vector.z:
			new_vector -= self.vector.project_z()
		elif vector.z < 0:
			new_vector += self.vector.project_z()

		# If vector is unchanged return
		# If vector has changed, constrain
		if new_vector == vector:
			return new_vector
		else:
			return self.constrain_vector_to_cell(new_vector)

	def set_basis_point(self, atomKey, position, interpolation='cubic'):
		"""Calculate and save basis function values for all points within cutoff region.

		Args:
			atomKey (int): Atom number, as given in Conquest_out
			position (Vector): 3D Cartesian real space vector
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
		"""
		atom = self.atoms[atomKey]
		Y = 0.0
		# Check if atom is in range
		relative_position = self.constrain_relative_vector(position - atom.atom_pos)
		if abs(relative_position) < self.atoms[atomKey].get_max_cutoff():
			for l in atom.radials:
				for zeta in atom.radials[l]:
					R = atom.get_radial_value_relative(l, zeta, relative_position, interpolation=interpolation)
					for m in range(-l, l+1):
						if R != 0:
							Y = sph(l, m, relative_position)
						self.basisPoints[atomKey][position][l][zeta][m] = R*Y
		else:
			self.basisPoints[atomKey][position] = SmartDict()

	def has_basis_point(self, atomKey, position):
		"""Check if basis point has already been calculated

		Args:
			atomKey (int): Atom number, as given in Conquest_out
			position (Vector): 3D Cartesian real space vector
		"""
		output = False
		if atomKey in self.basisPoints:
			if position in self.basisPoints[atomKey]:
				output = True
		return output

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
				self.bands[K].append(E)
			# Sort energy list
			self.bands[K] = sorted(self.bands[K])

	def get_gamma_energies(self):
		"""Return list of energies at gamma-point"""
		return sorted(self.bands[Vector.zero()])

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

	def get_psi(self, K, E, position, interpolation='cubic'):
		"""Evaluate wavefunction at specific position, k-point, and energy.

		Args:
			K (Vector): 3D Cartesian k-space vector
			E (float): Band energy
			position (Vector): 3D Cartesian real space vector
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'

		Returns:
			complex: Wavefunction value
		"""
		psi = complex(0.0, 0.0)
		# Loop over atoms
		for atomKey in self.atoms:
			# Check is basis function values have been calculated for this atom and point
			if not self.has_basis_point(atomKey, position):
				# Get and store basis function values
				# Store values so they only need to be calculated once
				self.set_basis_point(atomKey, position, interpolation=interpolation)
			# Use basis function value to calculate psi
			BP = self.basisPoints[atomKey][position]
			psi += self.atoms[atomKey].get_psi(K, E, position, basisPoint=BP, local=True)
		return psi

	def get_psi_gamma(self, E, position):
		"""Evaluate wavefunction at specific position and energy using only gamma-point.

		Args:
			E (float): Band energy
			position (Vector): 3D Cartesian real space vector

		Returns:
			complex: Wavefunction value
			"""
		return self.get_psi(Vector.zero(), E, position)

	def get_ldos(self, Emin, Emax, T, position, interpolation='cubic', debug=False):
		"""Evaluate local density of states (LDoS) within energy range at specific point.

		Args:
			Emin (float): Minimum energy
			Emax (float): Maximum energy
			T (float): Absolute temperature in K
			position (Vector): 3D Cartesian real space vector
			interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'
			debug (bool, opt.): Print extra information during runtime

		Returns:
			float: LDoS value
		"""
		I = 0.0
		totalK = len(self.bands)
		w = 1.0/totalK  # K-point weighting

		# Loop over all k-points
		for K in self.bands:
			for E in self.bands[K]:
				if Emin < E < Emax:
					# Calculate LDoS
					psi = self.get_psi(K, E, position, interpolation=interpolation)
					I += w * self.fermi_dirac(E, T) * (abs(psi)) ** 2
		if debug:
			print "Finished LDoS = "+str(I)+", at r = "+str(position)
		return I

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

	def calculate_support_grid(self, debug=False):
		"""Evaluate support function for each PAO on mesh.

		Args:
			debug (bool, opt.): If true, print extra information during runtime

		Returns:
			3D np array of SmartDict: Support function mesh, indexed by [x, y, z][atomKey][l][zeta][m]
		"""
		if debug:
			print "Calculating support grid"
		# Initialise support grid
		support_grid = np.empty_like(self.xMesh, dtype=SmartDict)

		# Iterate over all atoms
		for atomKey in self.atoms:
			atom = self.atoms[atomKey]

			# Get atom cutoff radius
			cut = atom.get_max_cutoff()

			# Get nearest mesh point to atom position
			atom_pos_on_mesh = self.get_nearest_mesh_vector(atom.atom_pos)

			# Get mesh points of maximum range of atoms orbitals in each direction
			x_lower_lim = self.get_nearest_mesh_value(atom.atom_pos.x - cut)
			x_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.x + cut)
			y_lower_lim = self.get_nearest_mesh_value(atom.atom_pos.y - cut)
			y_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.y + cut)
			z_lower_lim = self.get_nearest_mesh_value(atom.atom_pos.z - cut)
			z_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.z + cut)

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
						i = np.where(self.xMesh == constrained.x)[0][0]
						j = np.where(self.yMesh == constrained.y)[1][0]
						k = np.where(self.zMesh == constrained.z)[2][0]

						relative_position = self.constrain_relative_vector(r - atom_pos_on_mesh)
						# Iterate over orbitals
						for l in atom.radials:
							for zeta in atom.radials[l]:
								# Get radial part of wavefunction
								R = atom.get_radial_value_relative(l, zeta, relative_position)
								for m in range(-l, l + 1):
									# If R == 0, do not store
									if R != 0.0:
										# Get spherical harmonic
										Y = sph(l, m, relative_position)
										# Initialise support grid entry
										if not support_grid[i, j, k]:
											support_grid[i, j, k] = SmartDict()
										# Store support value
										support_grid[i, j, k][atomKey][l][zeta][m] = R * Y
			if debug:
				print "Calculated support grid for atom " + atom.ionName + " " + str(atomKey)
		return support_grid

	def calculate_ldos_grid(self, Emin, Emax, T, recalculate=False, write=True, debug=False):
		"""Calculate LDoS mesh.

		Args:
			Emin: Minimum energy
			Emax: Maximum energy
			T: Absolute temperature in Kelvin
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write calculated meshes to file
			debug (bool, opt.): Print extra information during runtime

		Returns:
			3D np.array: LDoS mesh
		"""
		if debug:
			print "Calculating LDoS grid"
		ldos_grid = np.zeros_like(self.xMesh, dtype=float)
		support_grid = self.get_support_grid(debug=debug, recalculate=recalculate)
		totalK = len(self.bands)
		w = 1.0 / totalK  # K-point weighting

		# Iterate over mesh points
		for i in range(self.xPoints):
			for j in range(self.yPoints):
				for k in range(self.zPoints):
					# Iterate over k-points
					if support_grid[i, j, k]:
						for K in self.bands:
							# Iterate over energies within range
							for E in self.bands[K]:
								if Emin <= E <= Emax:
									psi = complex(0.0, 0.0)
									# Iterate over atoms
									for atomKey in self.atoms:
										# Check if mesh point is with locality of atom
										if atomKey in support_grid[i, j, k]:
											atom = self.atoms[atomKey]
											# Iterate over orbitals
											for l in atom.bands[K][E]:
												for zeta in atom.bands[K][E][l]:
													for m in atom.bands[K][E][l][zeta]:
														# Calculate wavefunction
														coeff = atom.get_coefficient(K, E, l, zeta, m)
														psi += coeff * support_grid[i, j, k][atomKey][l][zeta][m]
									# Accumulate LDoS
									ldos_grid[i, j, k] += w * self.fermi_dirac(E, T) * (abs(psi)) ** 2
					if debug:
						print "Completed ", i, j, k
		return ldos_grid

	def get_support_grid(self, recalculate=False, write=True, debug=False):
		"""Get support function mesh.

		Args:
			recalculate (bool, opt.): Force recalculation, even if already stored
			write (bool, opt.): Write to file with name SUPPORT_FNAME
			debug (bool, opt.): Print extra information during runtime

		Returns:
			3D np.array of SmartDict: Support function mesh, indexed by [x, y, z][atomKey][l][zeta][m]
		"""
		# Read support grid from file if not stored by cell
		if not self.support_grid and not recalculate and os.path.isfile(SUPPORT_FNAME):
			self.support_grid = self.read_support_grid(debug=debug)
		else:
			# Recalculate support grid
			self.support_grid = self.calculate_support_grid(debug=debug)
			# Write to file
			if write:
				self.write_support(recalculate=False, debug=False)
		return self.support_grid

	def get_ldos_grid(self, Emin, Emax, T, recalculate=False, write=True, debug=False):
		"""Get LDoS mesh.

		Args:
			Emin: Minimum energy
			Emax: Maximum energy
			T: Absolute temperature in Kelvin
			recalculate (bool, opt.): Force recalculation of meshes, even if already stored
			write (bool, opt.): Write calculated grids to file
			debug (bool, opt.): Print extra information during runtime

		Returns:
			3D np.array: LDoS mesh
		"""
		# Read ldos grid from file if not stored by cell
		if not self.ldos_grid and not recalculate and os.path.isfile(LDOS_FNAME):
			self.ldos_grid = self.read_ldos_grid(debug=debug)
		else:
			# Calculate LDoS on mesh
			self.ldos_grid = self.calculate_ldos_grid(Emin, Emax, T, recalculate=recalculate, write=write, debug=debug)
			# Write to file
			if write:
				self.write_ldos(debug=debug, recalculate=False)
		return self.ldos_grid

	def write_support(self, recalculate=False, debug=False):
		"""Write support function mesh to file.

		Args:
			recalculate (bool, opt.): Force recalculation of meshes, even if already stored
			debug (bool, opt.): Print extra information during runtime
		"""
		# Get support mesh
		support_grid = self.get_support_grid(debug=debug, recalculate=recalculate)

		support_file = open(SUPPORT_FNAME, 'w')
		if debug:
			print "Writing support grid to "+SUPPORT_FNAME

		# Iterate over mesh points
		for i in range(self.xPoints):
			for j in range(self.yPoints):
				for k in range(self.zPoints):
					# If support function values exist at mesh point
					if support_grid[i, j, k]:
						support_file.write(str(i)+" "+str(j)+" "+str(k)+"\n")
						# Iterate over atoms
						for atomKey in support_grid[i, j, k]:
							# Write atom index
							support_file.write(str(atomKey)+"\n")
							# Iterate over orbitals
							for l in support_grid[i, j, k][atomKey]:
								for zeta in support_grid[i, j, k][atomKey][l]:
									for m in support_grid[i, j, k][atomKey][l][zeta]:
										# Write orbital data
										line = (str(l)+" "+str(zeta)+" "+str(m)+" "
												+str(support_grid[i, j, k][atomKey][l][zeta][m]))
										support_file.write(line+"\n")
		if debug:
			print "Support grid written to "+SUPPORT_FNAME
		support_file.close()

	def write_ldos(self, Emin, Emax, T, recalculate=False, debug=False):
		"""Write LDoS mesh to file.

		Args:
			Emin: Minimum energy
			Emax: Maximum energy
			T: Absolute temperature in K
			recalculate (bool, opt.): Force recalculation of mesh, even if stored
			debug (bool, opt.): Print extra information during runtime
		"""
		# Get LDoS mesh
		ldos_grid = self.get_ldos_grid(Emin, Emax, T, debug=debug, recalculate=recalculate)

		ldos_file = open(LDOS_FNAME, 'w')
		if debug:
			print "Write LDoS grid to "+LDOS_FNAME

		# Iterate over mesh points
		for i in range(self.xPoints):
			for j in range(self.yPoints):
				for k in range(self.zPoints):
					# If LDoS is non-zero at mesh point, write data to file
					if ldos_grid[i, j, k]:
						ldos_file.write(str(i)+" "+str(j)+" "+str(k)+" "+str(ldos_grid[i, j, k])+"\n")
		if debug:
			print "LDoS grid written to "+LDOS_FNAME
		ldos_file.close()

	def read_support_grid(self, debug=False):
		"""Read support function mesh from file"""
		support_file = open(SUPPORT_FNAME, 'r')
		support_grid = np.empty_like(self.xMesh, dtype=SmartDict)

		if debug:
			print "Reading support grid from "+SUPPORT_FNAME

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
						atomKey = int(line)
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
								support_grid[i, j, k][atomKey][l][zeta][m] = value
					else:
						reading_atoms = False
			except StopIteration:
				end_of_file = True
		if debug:
			print "Support grid successfully read"
		return support_grid

	def read_ldos_grid(self, debug=False):
		"""Read LDoS mesh from file"""
		ldos_file = open(LDOS_FNAME, 'r')
		ldos_grid = np.zeros_like(self.xMesh, dtype=float)

		if debug:
			print "Reading LDoS grid from "+LDOS_FNAME

		end_of_file = False
		line = ldos_file.next()
		line_split = line.split()
		while not end_of_file:
			try:
				# Get mesh indices
				i = int(line_split[0])
				j = int(line_split[1])
				k = int(line_split[2])
				# Get LDoS value
				value = float(line_split[3])
			except StopIteration:
				end_of_file = True
		if debug:
			print "LDoS grid successfully read"
		print np.max(ldos_grid)
		return ldos_grid
