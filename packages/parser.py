
import math
import os
import sys

import atomic
from cell import Cell
from smart_dict import SmartDict
from vector import Vector, KVector

__author__ = 'Johan G. McQuillan'
__email__ = 'johan.mcquillan.13@ucl.ac.uk'

HA_TO_EV = 27.2114  # Factor to convert Hartrees to electron volts


def get_ion(ion_name, ion_folder='ions'):
    """Create Ion object from .ion file."""
    
    try:
        # Open .ion
        ion_file = open(os.path.join(ion_folder, ion_name+'.ion'), 'r')

        # Skip preamble and first 9 lines
        line = ion_file.next()
        while '</preamble>' not in line:
            line = ion_file.next()
        for i in range(0, 9):
            ion_file.next()

        # Create empty Smartdict for Radial objects
        radial_dict = SmartDict()

        line = ion_file.next()
        while line.split()[0] != '#':
            # Read quantum numbers and zeta index
            metadata = line.split()
            l = int(metadata[0])
            n = int(metadata[1])
            zeta = int(metadata[2])

            # Get number of points for radial and cut-off radius
            line = ion_file.next()
            metadata = line.split()
            points = int(metadata[0])
            cutoff = float(metadata[2])

            # Initialise function data
            r = []
            R = []

            # Read data into Radial object and add to Ion
            for i in range(0, points):
                line = ion_file.next()
                x, y = line.split()
                x = float(x)
                # Multiply R by r^l
                y = float(y) * math.pow(x, l)
                r.append(x)
                R.append(y)

            # Create Radial object and store in SmartDict
            radial_dict[l][zeta] = atomic.Radial(n, l, zeta, r, R, cutoff)
            line = ion_file.next()
        ion_file.close()

        # Create Ion with radials
        return atomic.Ion(ion_name, radial_dict)

    except IOError:
        print '{} does not exist'.format(os.path.join(ion_folder, ion_name+'.ion'))
        sys.exit(1)


def get_cell(conquest_out, conquest_folder='conquest', ion_folder='ions', grid_spacing=0.5,
             group_size=400, weights=True, debug=False):
    """Return Cell object form CONQUEST files.

    Args:
        conquest_out (string): Base filename of simulation files.
            Within conquest_folder directory, there must exist files 'conquest_out',
            'conquest_out.dat', and 'conquest_out.dos'
        conquest_folder (string, opt): Folder for CONQUEST files.
        ion_folder (string, opt): Folder for .ion files.
        grid_spacing (float, opt): Mesh resolution in Bohr radii.
        group_size (int, opt): Maximum number of atoms to be saved to same support mesh file.
        weights (bool, opt): If true, look for k-point weightings in .dat file.
            Old version of Conquest_STMOutput did not provide k-point weights.
        debug (bool, opt): If true, print extra information during runtime.

    Returns:
        Cell: Simulation cell.
    """

    if debug:
        sys.stdout.write('Building simulation cell\n')
        sys.stdout.flush()
    try:
        # Open Conquest_out file
        conquest_out_file = open(os.path.join(conquest_folder, conquest_out), 'r')
        ions = {}
        atoms = {}

        # Skip opening lines
        line = conquest_out_file.next()
        while len(line.split()) != 8:
            line = conquest_out_file.next()

        atom_data = {}
        # Get atom positions and ion index
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
        while 'The simulation box has the following dimensions' not in line:
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

        # Get ion data
        while '------------------------------------------------------------------' not in line:
            raw_data = line.split()

            # Get ion index and name
            ion_type = int(raw_data[0])
            ion_name = raw_data[1]

            # Create Ion object
            ions[ion_name] = get_ion(ion_name, ion_folder=ion_folder)

            # Check for atoms of this ion type
            for atom_key in atom_data:
                atom_data_list = atom_data[atom_key]
                if atom_data_list[3] == ion_type:
                    # Get atom position
                    x, y, z = atom_data_list[:3]
                    r = Vector(x, y, z)

                    # Create Atom object
                    atoms[atom_key] = atomic.Atom(ion_name, r)

                    # Copy Ion attributes to Atom
                    atoms[atom_key].set_ion(ions[ion_name])

            conquest_out_file.next()
            line = conquest_out_file.next()

        conquest_out_file.close()

        try:
            # Open corresponding .dat file for basis coefficients
            conquest_dat_file = open(os.path.join(conquest_folder, conquest_out+'.dat'))
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

                    if weights:
                        # Get k-point weight
                        line = conquest_dat_file.next()
                        data = line.split()
                        K_weight = float(data[0])
                    else:
                        K_weight = 1

                    K = KVector(Kx, Ky, Kz, K_weight)

                    try:
                        line = conquest_dat_file.next()
                        while '#Kpoint' not in line:
                            # Get band energy
                            data = line.split()
                            bandN = int(data[0])
                            bandE = float(data[1]) * HA_TO_EV

                            # Get coefficient data
                            line = conquest_dat_file.next()
                            while len(line.split()) > 2:
                                data = line.split()
                                atom_index = int(data[0])
                                PAO = int(data[1])
                                coeff_string = data[2]

                                # Reformat string and create complex number
                                coeff_string = coeff_string.replace('(', '')
                                coeff_string = coeff_string.replace(')', '')
                                complex_string = coeff_string.split(',')
                                complex_coeff = complex(float(complex_string[0]),
                                                        float(complex_string[1]))

                                # Add coefficient to Atom
                                atoms[atom_index].add_coefficient(K, bandE, PAO, complex_coeff)
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
            print '{} does not exist'.format(os.path.join(conquest_folder, conquest_out+'.dat'))

        try:
            # Open .dos file
            conquest_dos_file = open(os.path.join(conquest_folder, conquest_out+'.dos'))
            line = conquest_dos_file.next()
            conquest_dos_file.close()

            # Get Fermi level
            data = line.split()
            fermi_lvl = float(data[2]) * HA_TO_EV

            # Create Cell
            cell = Cell(conquest_out, fermi_lvl, cell_length_x, cell_length_y, cell_length_z,
                        grid_spacing=grid_spacing, group_size=group_size)

            # Fill Cell with atoms
            for atom_key in atoms:
                cell.add_atom(atoms[atom_key], atom_key)

            return cell

        except IOError:
            print '{} does not exist'.format(os.path.join(conquest_folder, conquest_out+'.dos'))
            sys.exit(1)

    except IOError:
        print '{} does not exist'.format(os.path.join(conquest_folder, conquest_out))
        sys.exit(1)
