import math
import sys

import atomic
from cell import Cell
from smart_dict import SmartDict
from vector import Vector, KVector

HA_TO_EV = 27.2114  # Factor to convert Hartrees to electron volts


class Parser(object):
	"""Parses and stores data from input files"""

	def __init__(self, conquest_folder='conquest/', ion_folder='ions/'):
		"""Constructor for Parser"""
		self.conquest_folder = conquest_folder
		self.ion_folder = ion_folder

	def get_ion(self, ion_name):
		# Open .ion and initialise entry
		try:
			ion_file = open(self.ion_folder+ion_name+'.ion', 'r')

			# Skip preamble and first 9 lines
			line = ion_file.next()
			while '</preamble>' not in line:
				line = ion_file.next()
			for i in range(0, 9):
				ion_file.next()

			# Create empty dict of radial objects
			radial_dict = SmartDict()
			line = ion_file.next()
			while line.split()[0] != '#':
				# Read quantum numbers and zeta index
				metadata = line.split()
				l = int(metadata[0])
				n = int(metadata[1])
				zeta = int(metadata[2])

				# Get number of points for radial and cutoff radius
				line = ion_file.next()
				metadata = line.split()
				pts = int(metadata[0])
				cutoff = float(metadata[2])

				# Initialise R(r) function data
				r = []
				R = []
				# Read data into Radial object and add to Ion
				for i in range(0, pts):
					line = ion_file.next()
					x, y = line.split()
					x = float(x)
					y = float(y) * math.pow(x, l)
					r.append(x)
					R.append(y)

				# Create Radial object and store in dict
				radial_dict[l][zeta] = atomic.Radial(n, l, zeta, r, R, cutoff)
				line = ion_file.next()
			ion_file.close()
			# Create Ion with radials
			return atomic.Ion(ion_name, radial_dict)
		except IOError:
			print self.ion_folder+ion_name+'.ion does not exist'
			sys.exit(1)

	def make_cell(self, conquest_out, grid_spacing=0.5, group_size=150, debug=False):
		# Open Conquest_out file
		if debug:
			sys.stdout.write("Building simulation cell\n")
			sys.stdout.flush()
		try:
			conquest_out_file = open(self.conquest_folder+conquest_out, 'r')
			ions = {}
			atoms = {}

			# Skip opening lines
			line = conquest_out_file.next()
			while len(line.split()) != 8:
				line = conquest_out_file.next()

			# Read atomic positions
			atom_data = {}
			while len(line.split()) == 8 and line.split()[0].isdigit():
				raw_data = line.split()

				atom_index = int(raw_data[0])
				x = float(raw_data[1])
				y = float(raw_data[2])
				z = float(raw_data[3])
				ion_type = int(raw_data[4])

				atom_data[atom_index] = [x, y, z, ion_type]

				line = conquest_out_file.next()

			# Skip lines until cell data
			while "The simulation box has the following dimensions" not in line:
				line = conquest_out_file.next()
			line = conquest_out_file.next()

			# Get cell dimensions
			raw_data = line.split()
			cell_length_x = float(raw_data[2])
			cell_length_y = float(raw_data[5])
			cell_length_z = float(raw_data[8])

			# Skip lines until ion data
			while '------------------------------------------------------------------' not in line:
				line = conquest_out_file.next()
			conquest_out_file.next()
			conquest_out_file.next()
			line = conquest_out_file.next()

			# Get ion species
			while '------------------------------------------------------------------' not in line:
				raw_data = line.split()

				# Get ion data
				ion_type = int(raw_data[0])
				ion_name = raw_data[1]

				# Parse ion file and get Io object
				ions[ion_name] = self.get_ion(ion_name)

				# Check for atoms of this ion type
				for atom_key in atom_data:
					atom_data_list = atom_data[atom_key]
					if atom_data_list[3] == ion_type:
						# Get atom position
						x, y, z = atom_data_list[:3]
						r = Vector(x, y, z)

						# Create Atom object
						atoms[atom_key] = atomic.Atom(ion_name, r)

						# Apply Ion type to Atom
						atoms[atom_key].set_ion(ions[ion_name])
				conquest_out_file.next()
				line = conquest_out_file.next()
			conquest_out_file.close()

			# Open corresponding .dat file for basis coefficients
			try:
				conquest_dat_file = open(self.conquest_folder+conquest_out+'.dat')
				line = conquest_dat_file.next()

				# Loop over all lines
				end_of_file = False
				while not end_of_file:
					if '#Kpoint' in line:
						# Get k-point data
						line = conquest_dat_file.next()
						data = line.split()
						Kx = float(data[0])
						Ky = float(data[1])
						Kz = float(data[2])
						line = conquest_dat_file.next()
						data = line.split()
						K_weight = float(data[0])
						K = KVector(Kx, Ky, Kz, K_weight)
						try:
							line = conquest_dat_file.next()
							while '#Kpoint' not in line:
								# Get band energy data
								data = line.split()
								bandN = int(data[0])
								bandE = float(data[1]) * HA_TO_EV

								# Get coefficient data
								line = conquest_dat_file.next()
								while len(line.split()) > 2:
									data = line.split()
									a = int(data[0])
									PAO = int(data[1])
									coeff_string = data[2]
									coeff_string = coeff_string.replace('(', '')
									coeff_string = coeff_string.replace(')', '')
									complex_string = coeff_string.split(',')
									complex_coeff = complex(float(complex_string[0]), float(complex_string[1]))
									atoms[a].add_coefficient(K, bandE, PAO, complex_coeff)
									line = conquest_dat_file.next()
						except StopIteration:
							end_of_file = True
					else:
						# Check if end of file
						try:
							line = conquest_dat_file.next()
						except StopIteration:
							end_of_file = True
				conquest_dat_file.close()
			except IOError:
				print self.conquest_folder+conquest_out+'.dat does not exist'

			try:
				conquest_dos_file = open(self.conquest_folder+conquest_out+'.dos')
				line = conquest_dos_file.next()
				conquest_dos_file.close()

				data = line.split()
				fermi_lvl = float(data[2]) * HA_TO_EV

				cell = Cell(conquest_out, fermi_lvl, cell_length_x, cell_length_y, cell_length_z, grid_spacing=grid_spacing, group_size=group_size)

				for atom_key in atoms:
					atoms[atom_key].bands.lock()
					cell.add_atom(atoms[atom_key], atom_key)
				return cell
			except IOError:
				print self.conquest_folder+conquest_out+'.dos does not exist'
				sys.exit(1)
		except IOError:
			print self.conquest_folder+conquest_out+' does not exist'
			sys.exit(1)
