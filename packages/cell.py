import numpy as np
from copy import deepcopy

from sph import sph
from smartDict import SmartDict
from vector import Vector

BOLTZMANN = 8.6173303E-5  # Boltzmann's Constant in eV/K


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

	def __init__(self, name, fermiLevel, xLength, yLength, zLength, gridSpacing=0.5, debug=False):
		"""Constructs 3D cell with given dimensional.

		All lengths measured in Bohr radii (a0);
		All energies measured in Hartrees (Ha).

		Args:
			name (string): Name of simulation; used for plot titles
			fermiLevel (float): Fermi Level of simulation
			xLength (float): Length of cell along x
			yLength (float): Length of cell along y
			zLength (float): Length of cell along z
			gridSpacing (float, opt.): Resolution of mesh points
		"""

		self.name = name
		self.fermiLevel = fermiLevel
		self.gridSpacing = gridSpacing

		self.vector = Vector(xLength, yLength, zLength)

		# Calculator number of points
		self.xPoints = int(xLength / gridSpacing)
		self.yPoints = int(yLength / gridSpacing)
		self.zPoints = int(zLength / gridSpacing)

		# Form Cartesian meshes
		self.xMesh, self.yMesh, self.zMesh = np.mgrid[0: xLength: gridSpacing, 0: yLength: gridSpacing, 0: zLength: gridSpacing]

		# Initialise atoms and bands
		self.atoms = {}
		self.basisPoints = SmartDict()
		self.bands = {}
		self.support_grid = self.get_support_grid(debug=debug)

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

	def constrain_vector(self, vector):
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
			return self.constrain_vector(new_vector)

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
		relative_position = self.constrain_vector(position - self.atoms[atomKey].atom_pos)
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
			debug (bool, opt.): If true, print extra information during runtime

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

	def get_support_grid(self, debug=False):

		support_grid = SmartDict()
		for atomKey in self.atoms:
			atom = self.atoms[atomKey]
			support_grid[atomKey] = atom.get_grid(self.gridSpacing)
			if debug:
				print "Finished atom "+atom.ionName+" "+str(atomKey)
		return support_grid

	def
