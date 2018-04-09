
import os
import sys
import numpy as np
import warnings

from sph import sph
from smart_dict import SmartDict
from vector import Vector, KVector
from safe_open import safe_open

__author__ = 'Johan G. McQuillan'
__email__ = 'johan.mcquillan.13@ucl.ac.uk'


class Cell(object):
    """Simulation cell which holds Atom objects in a 3D mesh.

    All lengths measured in Bohr radii (a0).
    All energies measured in electron volts (eV).
    
    Attributes:
        name (string): Name of simulation; used for plot titles
        fermi_level (float): Fermi Level of simulation
        vector (Vector): Vector from origin to far corner of cell.
            Components are maximum lengths of cell in each dimension
        grid_spacing (float): Resolution of mesh points
        real_mesh (array(float)): Mesh of real space; [i, j, k] returns np.array([x, y, z])
        atoms ({int : Atom}): Atom objects of simulation indexed by atom number
        bands ({Vector : [float]}): List of band energies indexed by k-point vector
        support_mesh (array(SmartDict)): Mesh of support function values.
            Indexed by [i, j, k][atom_key][l][zeta][m]
    """

    BOLTZMANN = 8.6173303E-5        # Boltzmann's Constant in eV/K
    ELECTRON_MASS = 9.10938E-31     # Electron Mass in kg
    H_BAR = 4.135667662E-15         # Reduced Planck's Constant in eV.s
    DELTA_S_FACTOR = 2              # Default delta S given by minimum delta S scaled by this factor

    MESH_FOLDER = 'meshes/'         # Folder to save mesh files
    SUPPORT_FNAME = 'supp_'         # Prefix for support mesh files
    LDOS_FNAME = 'ldos_'            # Prefix for summed LDOS files
    PSI_FNAME = 'psi_'              # Prefix for wavefunction files
    PROP_PSI_FNAME = 'prop_'        # Prefix for propagated wavefunction files
    CURRENT_FNAME = 'current_'      # Prefix for current files
    EXT = '.dat'                    # Mesh file extension

    PRINT_RELATIVE_TO_EF = True     # Print energies as absolute or relative to Fermi level
    PROG_BAR_INTERVALS = 20         # Number of intervals in debug progress bar
    PROG_BAR_CHARACTER = '>'

    def __init__(self, name, fermi_level, x_length, y_length, z_length,
                 grid_spacing=0.5, group_size=400):
        """Constructs 3D cell with given dimensional.

        All lengths measured in Bohr radii (a0);
        All energies measured in Hartrees (Ha).

        Args:
            name (string): Name of simulation; used for plot titles.
            fermi_level (float): Fermi Level of simulation.
            x_length (float): Length of cell along x.
            y_length (float): Length of cell along y.
            z_length (float): Length of cell along z.
            grid_spacing (float, opt.): Resolution of mesh points.
            group_size (int, opt.): Maximum number of atoms to be saved to same support file.
        """

        self.name = name
        self.fermi_level = fermi_level
        self.grid_spacing = grid_spacing
        self.atom_group_size = group_size

        vector_x = int(x_length / grid_spacing) * grid_spacing
        vector_y = int(y_length / grid_spacing) * grid_spacing
        vector_z = int(z_length / grid_spacing) * grid_spacing
        self.vector = Vector(vector_x, vector_y, vector_z)

        # Initialise Cartesian meshes
        self.real_mesh = np.transpose(np.mgrid[0:x_length:grid_spacing,
                                               0: y_length: grid_spacing, 0:
                                               z_length: grid_spacing],
                                      (1, 2, 3, 0))

        # Total number of mesh points
        self.mesh_points = self.real_mesh.size

        # Initialise atoms and bands
        self.atoms = {}
        self.bands = {}

        # Currently support function mesh
        self.support_mesh = None
        # Support meshes are split into groups, each with self.group_size atoms
        # For number of atoms > 400, storing mesh for all atoms takes huge amounts of RAM
        # First group is 0, second is 1 etc.
        # -1 indicates no group stored
        self.current_group = -1  # Current atom group stored

        self.default_delta_s = self.delta_s()  # Default delta_s

        # Vectorised method to calculate wavefunction
        self.psi_vec = np.vectorize(self.calculate_psi_grid_vec)

    def energy_list(self):
        """Return sorted list of energies from all k-points."""
        
        energies = []
        for K in self.bands:
            energies.extend(self.bands[K])
        return sorted(energies)

    def has_band(self, K, E):
        """Check if cell stores specified band.

        Args:
            K (Vector): 3D Cartesian k-space vector.
            E (float): Band energy.
        """
        
        output = False
        if K in self.bands:
            if E in self.bands[K]:
                output = True
        return output

    def add_atom(self, atom, atom_key):
        """Add atom to self.atoms, indexed by atom_key.

        Args:
            atom (Atom): Atom object.
            atom_key (int): Atom number, as given in Conquest_out.
        """
        
        # Add atom to dict
        self.atoms[atom_key] = atom

        # Add band energies and k-points to self.bands
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
            energy (float): Energy in eV.
            temperature (float): Absolute temperature in K.

        Returns:
            float: Occupation factor, between 0 and 1.
        """
        
        with warnings.catch_warnings():
            # Suppress RuntimeWarning from overflow and underflow in np.exp
            warnings.simplefilter('ignore', RuntimeWarning)
            f = 1.0 / (np.exp((energy - self.fermi_level) / (self.BOLTZMANN * temperature)) + 1)
        return f

    def bias_to_energy_range(self, V):
        """Return absolute energy range from the Fermi level to the bias voltage."""
        
        if V > 0:
            min_E = self.fermi_level
            max_E = self.fermi_level + V
        else:
            min_E = self.fermi_level + V
            max_E = self.fermi_level
        return min_E, max_E

    def get_nearest_mesh_value(self, x, points=None):
        """Return nearest mesh point to x. Not constrained to lie within simulation cell.

        Args:
            x (float): Value in a0 to find nearest point; Works for x, y, and z dimensions.
            points (float, opt.): Number of points in cell dimension.
                If given, will return mesh index of the found point

        Returns:
            float: Nearest mesh point value
            int: Index of mesh point; only given if points argument is specified
        """
        
        # Get quotient and remainder wrt grid spacing
        div_x, mod_x = divmod(x, self.grid_spacing)
        # Check if x should be rounded up or down
        if mod_x >= self.grid_spacing / 2:
            div_x += 1
        # Get new point
        new_x = div_x * self.grid_spacing
        if points is not None:
            i = div_x - 1
            while i < 0:
                i += points
            while i >= points:
                i -= points
            return new_x, int(i)
        else:
            return new_x

    def constrain_relative_vector(self, vector):
        """Using periodic boundaries, return smallest Vector that is equivalent to input Vector."""
        
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

    def calculate_support_group(self, group, interpolation='cubic', debug=False):
        """Evaluate support function for each PAO for a given atom group.

        Args:
            group (int): Atom group
            interpolation (string, opt.): Method of interpolation.
                Possible arguments are 'linear', 'quadratic', 'cubic'.
            debug (bool, opt.): If true, print extra information during runtime.

        Returns:
            array(SmartDict): Mesh of support function values.
                Indexed by [x, y, z][atom_key][l][zeta][m]
        """
        
        if debug:
            print 'Calculating support mesh for atom group', str(group)
        # Initialise support grid
        support_grid = np.empty(self.real_mesh.shape[:3], dtype=SmartDict)

        # Get first and last atom numbers for this group
        lower_bound = group * self.atom_group_size
        upper_bound = (group + 1) * self.atom_group_size

        # Get atom keys for group
        atom_list = sorted([a for a in self.atoms.keys() if lower_bound <= a < upper_bound])

        # Iterate over atoms
        for atom_key in atom_list:
            atom = self.atoms[atom_key]

            debug_str = '  Support grid for atom {:3} - {!s}: '.format(atom_key, atom.ion_name)
            if debug:
                sys.stdout.write(debug_str)
                sys.stdout.flush()

            # Get atom cutoff radius
            cut = atom.get_max_cutoff()

            # Get mesh points of maximum range of atoms orbitals in each direction
            x_lower_lim, i_start = self.get_nearest_mesh_value(atom.atom_pos.x - cut,
                                                               points=self.real_mesh.shape[0])
            x_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.x + cut) + self.grid_spacing
            y_lower_lim, j_start = self.get_nearest_mesh_value(atom.atom_pos.y - cut,
                                                               points=self.real_mesh.shape[1])
            y_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.y + cut) + self.grid_spacing
            z_lower_lim, k_start = self.get_nearest_mesh_value(atom.atom_pos.z - cut,
                                                               points=self.real_mesh.shape[2])
            z_upper_lim = self.get_nearest_mesh_value(atom.atom_pos.z + cut) + self.grid_spacing

            # Get array of mesh points within cutoff
            local_mesh = np.transpose(np.mgrid[x_lower_lim:x_upper_lim:self.grid_spacing,
                                               y_lower_lim:y_upper_lim:self.grid_spacing,
                                               z_lower_lim:z_upper_lim:self.grid_spacing],
                                      (1, 2, 3, 0))
            lm_shape = local_mesh.shape[:3]

            # Progress bar initialisation
            points_done = 0
            bars_done = 0
            total_points = local_mesh.shape[0]*local_mesh.shape[1]*local_mesh.shape[2]

            # The local mesh may exist over the periodic boundaries
            # Roll the full support mesh such the [0, 0, 0] entry corresponds
            #   physically to the same point as the [0, 0, 0] point on the local mesh
            rolled_mesh = np.roll(support_grid, -i_start, 0)
            rolled_mesh = np.roll(rolled_mesh, -j_start, 1)
            rolled_mesh = np.roll(rolled_mesh, -k_start, 2)

            # Extract the part of the full support mesh that corresponds to
            #   the same physical space as the local mesh
            partial_mesh = rolled_mesh[0:lm_shape[0], 0:lm_shape[1], 0:lm_shape[2]]

            # Iterate over the local mesh
            for local_ijk in np.ndindex(lm_shape):
                position = local_mesh[local_ijk]
                r = Vector(*position)

                # Find the shortest Vector between r and atom_pos using periodic boundaries
                relative_position = self.constrain_relative_vector(r - atom.atom_pos)

                # Convert indices of local mesh point from tuple into a list
                partial_ijk_list = list(local_ijk)
                for i in range(len(partial_ijk_list)):
                    # If point is outside boundary, reduce to lie within cell
                    if partial_ijk_list[i] >= self.real_mesh.shape[i]:
                        partial_ijk_list[i] -= self.real_mesh.shape[i]
                # Convert back to tuple
                partial_ijk = tuple(partial_ijk_list)

                # Iterate over orbitals
                for l in atom.radials:
                    for zeta in atom.radials[l]:
                        # Get radial part of wavefunction
                        R = atom.get_radial_value_relative(l, zeta, relative_position,
                                                           interpolation=interpolation)

                        # If R == 0, do not store
                        if R != 0.0:
                            for m in range(-l, l + 1):
                                # Get spherical harmonic
                                Y = sph(l, m, relative_position)

                                # Initialise support mesh entry
                                if partial_mesh[partial_ijk] is None:
                                    partial_mesh[partial_ijk] = SmartDict()

                                if m not in partial_mesh[partial_ijk][atom_key][l][zeta]:
                                    partial_mesh[partial_ijk][atom_key][l][zeta][m] = 0.0

                                # Store support function value
                                partial_mesh[partial_ijk][atom_key][l][zeta][m] += R * Y
                points_done += 1

                # Update progress bar
                prog = float(points_done) / total_points
                if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
                    percent = prog * 100
                    sys.stdout.write('\r')
                    sys.stdout.write(debug_str)
                    sys.stdout.write(
                        ' [{:<{}}]'.format(self.PROG_BAR_CHARACTER * bars_done,
                                           self.PROG_BAR_INTERVALS))
                    sys.stdout.write(' {:3.0f}%'.format(percent))
                    sys.stdout.flush()
                    bars_done += 1

            # Copy partial support mesh into full support mesh
            rolled_mesh[0:lm_shape[0], 0:lm_shape[1], 0:lm_shape[2]] = partial_mesh

            # Roll mesh back to original position
            rolled_mesh = np.roll(rolled_mesh, i_start, 0)
            rolled_mesh = np.roll(rolled_mesh, j_start, 1)
            support_grid = np.roll(rolled_mesh, k_start, 2)

            if debug:
                sys.stdout.write('\n')
                sys.stdout.flush()

        return support_grid

    def support_group_filename(self, group):
        """Return standardised filename for relevant support function file."""
        
        return self.MESH_FOLDER+self.SUPPORT_FNAME+self.name+'_'+str(self.grid_spacing)+'_'+str(self.atom_group_size)+'_'+str(group)+self.EXT

    def write_support_group(self, group, support_mesh, debug=False):
        """Write support function group to file."""
        
        filename = self.support_group_filename(group)
        support_file = safe_open(filename, 'w')

        if debug:
            sys.stdout.write('Writing support group to '+filename+': ')
            sys.stdout.flush()
        points_done = 0
        bars_done = 0

        # Iterate over mesh points
        for ijk in np.ndindex(self.real_mesh.shape[:3]):
            i, j, k = ijk

            # If support function values exist at mesh point
            if support_mesh[ijk]:
                support_file.write(str(i) + ' ' + str(j) + ' ' + str(k) + '\n')

                # Iterate over atoms
                for atom_key in support_mesh[ijk]:
                    # Write atom index
                    support_file.write(str(atom_key) + '\n')

                    # Iterate over orbitals
                    for l in support_mesh[ijk][atom_key]:
                        for zeta in support_mesh[ijk][atom_key][l]:
                            for m in support_mesh[ijk][atom_key][l][zeta]:
                                # Write orbital data
                                line = (str(l) + ' ' + str(zeta) + ' ' + str(m) + ' '
                                        + str(support_mesh[ijk][atom_key][l][zeta][m]))
                                support_file.write(line + '\n')

            # Update progress bar
            points_done += 1
            if debug and float(points_done) / self.mesh_points * self.PROG_BAR_INTERVALS > bars_done:
                sys.stdout.write('\r')
                sys.stdout.write(' [{:<{}}]'.format(self.PROG_BAR_CHARACTER * bars_done,
                                                    self.PROG_BAR_INTERVALS))
                sys.stdout.flush()
                bars_done += 1

        if debug:
            sys.stdout.write('\n')
            sys.stdout.flush()
        support_file.close()

    def read_support_group(self, group, debug=False):
        """Read support function group from file."""
        
        filename = self.support_group_filename(group)
        support_file = open(filename, 'r')
        support_mesh = np.empty(self.real_mesh.shape[:3], dtype=SmartDict)

        if debug:
            print 'Reading support group from '+filename

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
            print 'Support grid successfully read'
        return support_mesh

    def get_support_group(self, group, recalculate=False, interpolation='cubic', debug=False):
        """Get support function mesh for given atom group.

        If it is already stored in a file and recalculate is False, then it will read from file.
        Otherwise, it will calculate support functions.

        Args:
            group (int): Atom group
            recalculate (bool, opt.): Force recalculation, even if already stored.
            interpolation (string, opt.): Method of interpolation.
                Possible arguments are 'linear', 'quadratic', 'cubic'.
            debug (bool, opt.): Print extra information during runtime.

        Returns:
            array(SmartDict): Mesh of support function values.
                Indexed by [i, j, k][atom_key][l][zeta][m].
        """
        
        # Check if support function group is already stored in field
        if group == self.current_group and self.support_mesh is not None:
            return self.support_mesh
        else:
            if not recalculate and os.path.isfile(self.support_group_filename(group)):
                # Read support grid from file
                support_mesh = self.read_support_group(group, debug=debug)
            else:
                # Recalculate support grid
                support_mesh = self.calculate_support_group(group, interpolation=interpolation,
                                                            debug=debug)
                # Write to file
                self.write_support_group(group, support_mesh, debug=debug)
            self.current_group = group
            self.support_mesh = support_mesh
            return support_mesh

    def calculate_psi_grid_vec(self, support_dict, atom_key, l, zeta, m, coefficient):
        """Calculate wavefunction contribution from a point of the support function mesh
        for a given orbital and atom."""
        
        if (support_dict is not None and atom_key in support_dict and l in support_dict[atom_key]
                and zeta in support_dict[atom_key][l] and m in support_dict[atom_key][l][zeta]):
            psi = coefficient * support_dict[atom_key][l][zeta][m]
        else:
            psi = complex(0, 0)
        return psi

    def calculate_psi_grid(self, K, E, recalculate=False, vectorised=True, interpolation='cubic',
                           debug=False):
        """Evaluate wavefunction over the mesh for a given k-point and energy.

        Args:
            K (KVector): K-point Vector.
            E (float): Band energy.
            recalculate (bool, opt.): Force recalculation, even if already stored.
            vectorised (bool, opt.): If true, use NumPy vectorisation.
            interpolation (string, opt.): Method of interpolation.
                Possible arguments are 'linear', 'quadratic', 'cubic'
            debug (bool, opt.): Print extra information during runtime.

        Returns:
            array(complex): Mesh of complex wavefunction values.
        """

        # Initialise progress bar
        atoms_done = 0
        total_atoms = len(self.atoms)
        bars_done = 0

        # Initialise wavefunction mesh
        psi_grid = np.zeros_like(self.real_mesh[..., 0], dtype=complex)

        # Print debug info
        if self.PRINT_RELATIVE_TO_EF:
            E_str = str(E - self.fermi_level) + ' eV'
        else:
            E_str = str(E) + ' eV'
        debug_str = 'Calculating psi(r) at k = '+str(K)+', E = '+E_str+': '
        if debug:
            sys.stdout.write(debug_str)
            sys.stdout.flush()

        # Iterate over all atoms
        previous_group = 0
        support_grid = self.get_support_group(0, recalculate=recalculate, debug=debug)
        for atom_key in sorted(self.atoms.iterkeys()):
            atom = self.atoms[atom_key]
            group = atom_key / self.atom_group_size

            # Check if current atom group is the one currently read from file
            if group != previous_group:
                support_grid = self.get_support_group(group, recalculate=recalculate,
                                                      interpolation=interpolation, debug=debug)
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
            # Update progress bar
            atoms_done += 1
            prog = float(atoms_done) / total_atoms
            if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
                percent = prog * 100
                sys.stdout.write('\r')
                sys.stdout.write(debug_str)
                sys.stdout.write(' [{:<{}}]'.format(self.PROG_BAR_CHARACTER * bars_done,
                                                    self.PROG_BAR_INTERVALS))
                sys.stdout.write(' {:3.0f}%'.format(percent))
                sys.stdout.flush()
                bars_done = int(prog * self.PROG_BAR_INTERVALS)
        if debug:
            sys.stdout.write('\n')
            sys.stdout.flush()
        return psi_grid

    def psi_filename(self, K, E):
        """Return standardised filename for relevant wavefunction file."""
        
        return (self.MESH_FOLDER+self.PSI_FNAME+self.name+'_'+str(self.grid_spacing)+'_'
                +str(K.x)+'_'+str(K.y)+'_'+str(K.z)+'_'+str(E)+self.EXT)

    def write_psi_grid(self, psi_grid, K, E):
        """Write wavefunction function mesh to file."""
        
        filename = self.psi_filename(K, E)
        psi_file = safe_open(filename, 'w')
        for ijk in np.ndindex(self.real_mesh.shape[:3]):
            i, j, k = ijk
            if psi_grid[ijk] != 0:
                psi = psi_grid[ijk]
                psi_file.write(str(i)+' '+str(j)+' '+str(k)+' '+str(psi.real)+' '+str(psi.imag)+'\n')
        psi_file.close()

    def read_psi_grid(self, K, E, debug=False):
        """Read wavefunction mesh from file."""
        
        filename = self.psi_filename(K, E)
        psi_file = open(filename, 'r')
        psi_grid = np.zeros_like(self.real_mesh[..., 0], dtype=complex)

        if debug:
            if self.PRINT_RELATIVE_TO_EF:
                E_str = str(E - self.fermi_level) + ' eV'
            else:
                E_str = str(E) + ' eV'
            sys.stdout.write('Reading psi(r) at k = {!s}, E = {}\n'.format(K, E_str))

        for line in psi_file:
            line_split = line.split()
            i = int(line_split[0])
            j = int(line_split[1])
            k = int(line_split[2])
            real = float(line_split[3])
            imag = float(line_split[4])
            psi_grid[i, j, k] = complex(real, imag)
        return psi_grid

    def get_psi_grid(self, K, E, recalculate=False, vectorised=True, interpolation='cubic',
                     debug=False):
        """Get mesh of complex wavefunction values.

        If it is already stored in a file and recalculate is False, then it will read from file.
        Otherwise, it will calculate wavefunction.

        Args:
            K (KVector): K-point Vector
            E (float): Band energy
            recalculate (bool, opt.): Force recalculation, even if already stored
            vectorised (bool, opt.): If true, use NumPy vectorisation
            interpolation (string, opt.): Method of interpolation.
                Possible arguments are 'linear', 'quadratic', 'cubic'
            debug (bool, opt.): Print extra information during runtime

        Returns:
            array(complex): Mesh of complex wavefunction values
        """
        
        if not recalculate and os.path.isfile(self.psi_filename(K, E)):
            # Read data from file
            psi_grid = self.read_psi_grid(K, E, debug=debug)
        else:
            psi_grid = self.calculate_psi_grid(K, E, recalculate=recalculate, vectorised=vectorised,
                                               interpolation=interpolation, debug=debug)
            self.write_psi_grid(psi_grid, K, E)
        return psi_grid

    def calculate_ldos_grid(self, min_E, max_E, T, recalculate=False, vectorised=True,
                            interpolation='cubic', debug=False):
        """Calculate summed LDOS.

        Args:
            min_E (float): Minimum absolute energy
            max_E (float): Maximum absolute energy
            T (float): Absolute temperature in Kelvin
            recalculate (bool, opt.): Force recalculation, even if already stored
            vectorised (bool, opt.): If true, use NumPy vectorisation
            interpolation (string, opt.): Method of interpolation.
                Possible arguments are 'linear', 'quadratic', 'cubic'
            debug (bool, opt.): Print extra information during runtime

        Returns:
            array(float): LDoS mesh
        """
        
        # Print debug info
        if debug:
            sys.stdout.write('Calculating local density of states grid\n')
            sys.stdout.flush()

        # Initialise mesh
        ldos_grid = np.zeros_like(self.real_mesh[..., 0], dtype=float)

        total_k_weight = 0
        for K in self.bands:
            total_k_weight += K.weight

        # Iterate over energies
        for K in self.bands:
            for E in self.bands[K]:
                if min_E <= E <= max_E:
                    psi_grid = self.get_psi_grid(K, E, recalculate=recalculate,
                                                 vectorised=vectorised, interpolation=interpolation,
                                                 debug=debug)
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
        
        return self.MESH_FOLDER+self.LDOS_FNAME+self.name+'_'+str(self.grid_spacing)+'_'+str(min_E)+\
            '_'+str(max_E)+'_'+str(T)+self.EXT

    def write_ldos_grid(self, ldos_grid, min_E, max_E, T, debug=False):
        """Write LDoS mesh to file.

        Args:
            ldos_grid (array(float)): 3D array of LDOS values.
            min_E (float): Minimum energy.
            max_E (float): Maximum energy.
            T (float): Absolute temperature in K.
            debug (bool, opt.): Print extra information during runtime.
        """
        
        filename = self.ldos_filename(min_E, max_E, T)
        # Get LDOS mesh
        ldos_file = safe_open(filename, 'w')
        if debug:
            sys.stdout.write('Writing LDOS grid to {}\n'.format(filename))
            sys.stdout.flush()
        # Iterate over mesh points
        for i, j, k in np.ndindex(self.real_mesh.shape[:3]):
            # If LDOS is non-zero at mesh point, write data to file
            if ldos_grid[i, j, k]:
                ldos_file.write(str(i)+' '+str(j)+' '+str(k)+' '+str(ldos_grid[i, j, k])+'\n')
        ldos_file.close()

    def read_ldos_grid(self, min_E, max_E, T, debug=False):
        """Read LDOS mesh from file."""
        
        filename = self.ldos_filename(min_E, max_E, T)
        ldos_file = open(filename, 'r')
        ldos_grid = np.zeros_like(self.real_mesh[..., 0])

        if debug:
            print 'Reading LDoS grid from '+filename

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
            print 'LDoS grid successfully read'
        return ldos_grid

    def get_ldos_grid(self, min_E, max_E, T, recalculate=False, vectorised=True,
                      interpolation='cubic', debug=False):
        """Get LDOS mesh by calculating or reading from file.

        Args:
            min_E (float): Minimum absolute energy
            max_E (float): Maximum absolute energy
            T (float): Absolute temperature in K
            recalculate (bool, opt.): Force recalculation of meshes, even if already stored
            vectorised (bool, opt.): If true, use NumPy vectorisation
            interpolation (string, opt.): Method of interpolation.
                Possible arguments are 'linear', 'quadratic', 'cubic'
            debug (bool, opt.): Print extra information during runtime

        Returns:
            array(float): 3D mesh of LDOS values
        """
        
        # Read ldos grid from file if not stored by cell
        if not recalculate and os.path.isfile(self.ldos_filename(min_E, max_E, T)):
            ldos_grid = self.read_ldos_grid(min_E, max_E, T, debug=debug)
        else:
            # Calculate LDOS on mesh
            ldos_grid = self.calculate_ldos_grid(min_E, max_E, T, recalculate=recalculate,
                                                 vectorised=vectorised, interpolation=interpolation,
                                                 debug=debug)
            self.write_ldos_grid(ldos_grid, min_E, max_E, T, debug=debug)
        return ldos_grid

    def periodic_gradient(self, mesh):
        """Calculate gradient of mesh with periodic boundary conditions enforced.

        This is done by padding the mesh with extra columns in each dimension, then copying the
        columns from the other side of the mesh into the pads. This makes the columns corresponding
        to a boundary of the simulation cell effectively next to the columns near the opposing
        boundary, such that numerical calculation of the gradient should be continuous over the
        boundary.
        """

        # Width of padding
        pad = 3
        # Get shape of padded array
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
        padded_gradient = np.transpose(np.array(np.gradient(padded_mesh, self.grid_spacing)),
                                       (1, 2, 3, 0))

        # Return unpadded gradient
        return padded_gradient[pad:-pad, pad:-pad, pad:-pad]

    def delta_s(self):
        """Calculate the default value for delta S.
        
        This is the minimum delta S from Paz and Soler, multiplied by DELTA_S_FACTOR."""
        
        min_delta_s = 2.0 * self.grid_spacing / self.H_BAR * np.sqrt(4.85*2.0*self.ELECTRON_MASS)
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

    def broadened_surface(self, charge_density_mesh, fraction, max_height_index, delta_s=None):
        """Calculate the magnitude of the c mesh"""
        if delta_s is None:
            delta_s = self.default_delta_s

        # Get isovalue
        max_value = np.max(charge_density_mesh)
        isovalue = fraction * max_value

        # Evaluate logarithm on all non-zero entrie
        log_mesh = np.empty(charge_density_mesh.shape, dtype=float)
        zero_points = charge_density_mesh == 0
        log_mesh[~zero_points] = np.log(charge_density_mesh[charge_density_mesh != 0] / isovalue)
        log_mesh[zero_points] = np.inf

        # Apply broadening to surface
        broadened_mesh = np.where(abs(log_mesh) < delta_s,
                                  15.0 / (16.0 * delta_s) * (1.0 - (log_mesh / delta_s) ** 2) ** 2, 0)

        # Remove overlapping layers of surface
        # Iterate over x and y
        for i, j in np.ndindex(broadened_mesh.shape[:2]):
            on_surface = False
            past_surface = False
            # Iterate over z, starting from top
            for k in reversed(range(broadened_mesh.shape[2])):
                if k < max_height_index:
                    if past_surface:
                        # First surface has been traversed, replace subsequent elements with zeros
                        broadened_mesh[i, j, k] = 0
                    elif broadened_mesh[i, j, k] != 0 and not on_surface:
                        # Highest surface has been reached
                        on_surface = True
                    elif on_surface and broadened_mesh[i, j, k] == 0:
                        # Was on surface, now just below
                        on_surface = False
                        past_surface = True
                else:
                    broadened_mesh[i, j, k] = 0
        return broadened_mesh

    def get_c(self, charge_density_mesh, fraction, tip_height_index, delta_s=None):
        """Return c mesh for surface integration.

        Args:
            charge_density_mesh (array(float)): Charge density or LDOS mesh.
            fraction (float): Isovalue given by fraction of maximum mesh value.
            tip_height_index (int): Z-index of tip height.
            delta_s (float, opt.): Durface broadening parameter; If None, uses default value.
        """
        
        if delta_s is None:
            delta_s = self.default_delta_s

        # Find all points of zero density
        zero_points = charge_density_mesh == 0

        # Get magnitude of c mesh
        broadened_mesh = self.broadened_surface(charge_density_mesh, fraction, tip_height_index,
                                                delta_s=delta_s)

        # Get gradient of density mesh
        gradient_surface = self.periodic_gradient(charge_density_mesh)

        # Calculate final mesh
        vector_surface = np.zeros(gradient_surface.shape, dtype=float)
        for i in range(3):
            gradient_surface[~zero_points, i] = gradient_surface[~zero_points, i] / charge_density_mesh[~zero_points]
            gradient_surface[zero_points, i] = 0
            vector_surface[..., i] = np.multiply(broadened_mesh, gradient_surface[..., i])

        return vector_surface

    def get_A_mesh(self, c, wavefunction_mesh):
        """Calculate A mesh."""
        
        grad_wavefunction = self.periodic_gradient(wavefunction_mesh)
        return self.mesh_dot_product(c, grad_wavefunction)

    @staticmethod
    def mesh_dot_product(vector_mesh_A, vector_mesh_B):
        """Return dot/scalar product of two vector meshes."""
        
        # Check if resulting mesh will be complex
        if vector_mesh_A.dtype == complex or vector_mesh_B.dtype == complex:
            dtype = complex
        else:
            dtype = float

        scalar_mesh = np.zeros(vector_mesh_A.shape[:3], dtype=dtype)
        for i in range(3):
            scalar_mesh += vector_mesh_A[..., i]*vector_mesh_B[..., i]
        return scalar_mesh

    def get_B_mesh(self, c, wavefunction_mesh):
        """Calculate B mesh."""
        
        B = np.zeros_like(c, dtype=complex)
        for i in range(3):
            B[..., i] = c[..., i] * wavefunction_mesh
        return B

    def greens_function_mesh(self, z_index, tip_work_func, tip_energy, debug=False):
        """Calculate Tersoff-Hamann tip Greens function over entire mesh.

        Uses symmetry of the function to reduce number of calculations.

        Args:
            z_index (float): z-index of tip plane.
            tip_work_func (float): Work function of tip.
            tip_energy (float): Fermi-level of tip.
            debug (bool, opt.): Print extra information during runtime.
        """

        # Create square grid with sides equal to the size of the x or y=axis of the real mesh
        #   whichever is larger
        plane_length_half = - (- (max(self.real_mesh.shape[:2]) + 1) / 2)
        plane_length_full = plane_length_half * 2 - 1
        plane_shape_half = (plane_length_half,) * 2
        plane_UR = np.zeros(plane_shape_half, dtype=float)

        # Initialise Green's function mesh
        G_shape = (plane_length_full, plane_length_full, self.real_mesh.shape[2])
        G_mesh = np.zeros(G_shape, dtype=float)

        # Print debug info
        debug_str = 'Calculating G(r - R): '
        if debug:
            sys.stdout.write(debug_str)
            sys.stdout.flush()
        points_done = 0
        bars_done = 0

        # Iterate over z values
        for k in range(G_shape[2]):
            # Iterate over right-angle triangle in x and y
            # The tip is considered to be in corner of the mesh at i = j = 0
            for i in range(plane_length_half):
                for j in range(i + 1):
                    # If above tip, define as 0
                    if k >= z_index:
                        G = 0
                    else:
                        # Calculate Green's function
                        distance = self.grid_spacing * np.sqrt(i**2 + j**2 + (z_index - k)**2)
                        G = self.greens_function(distance, tip_work_func, tip_energy)

                    # Copy values from this side of triangle into other
                    # In effect, reflect values along x = y
                    for x, y in [i, j], [j, i]:
                        plane_UR[x, y] = G
            # Flip square to create four subsquares
            plane_UL = np.flipud(plane_UR[1:, :])
            plane_LR = np.fliplr(plane_UR[:, 1:])
            plane_LL = np.flipud(plane_LR[1:, :])

            # Combine subsquares
            # Tip is now defined above the centre of mesh
            plane_U = np.concatenate((plane_UL, plane_UR), axis=0)
            plane_L = np.concatenate((plane_LL, plane_LR), axis=0)
            plane = np.concatenate((plane_L, plane_U), axis=1)

            G_mesh[..., k] = plane

            # Update progress bar
            points_done += 1
            prog = float(points_done) / self.real_mesh.shape[2]
            if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
                percent = prog * 100
                sys.stdout.write('\r')
                sys.stdout.write(debug_str)
                sys.stdout.write(' [{:<{}}]'.format(self.PROG_BAR_CHARACTER * bars_done,
                                                    self.PROG_BAR_INTERVALS))
                sys.stdout.write(' {:3.0f}%'.format(percent))
                sys.stdout.flush()
                bars_done += 1

        if debug:
            sys.stdout.write('\n')
            sys.stdout.flush()
        return G_mesh

    def propagated_psi_filename(self, K, E, T, fraction, z, delta_s=None):
        """Return standardised filename for relevant propagated wavefunction file."""
        
        if delta_s is None:
            delta_s = self.default_delta_s
        return (self.MESH_FOLDER+self.PROP_PSI_FNAME+self.name+'_'+str(self.grid_spacing)+'_'
                +str(K.x)+'_'+str(K.y)+'_'+str(K.z)+'_'+str(E)+'_'+str(T)+'_'+str(fraction)+'_'
                +str(z)+'_'+str(delta_s)+self.EXT)

    def write_prop_psi(self, psi, K, E, T, fraction, z, delta_s=None):
        """Write propagated wavefunction function mesh to file."""
        
        if delta_s is None:
            delta_s = self.default_delta_s
        filename = self.propagated_psi_filename(K, E, T, fraction, z, delta_s=delta_s)
        psi_file = safe_open(filename, 'w')
        for i, j in np.ndindex(psi.shape):
            if psi[i, j] != 0:
                p = psi[i, j]
                psi_file.write(str(i)+' '+str(j)+' '+str(p.real)+' '+str(p.imag)+'\n')
        psi_file.close()

    def read_prop_psi(self, K, E, T, fraction, z, delta_s=None, debug=False):
        """Read propagated wavefunction mesh from file."""
        
        if delta_s is None:
            delta_s = self.default_delta_s

        filename = self.propagated_psi_filename(K, E, T, fraction, z, delta_s=delta_s)
        psi_file = open(filename, 'r')
        psi = np.zeros(self.real_mesh.shape[:2], dtype=complex)

        if debug:
            if self.PRINT_RELATIVE_TO_EF:
                E_str = str(E - self.fermi_level) + ' eV'
            else:
                E_str = str(E) + ' eV'
            sys.stdout.write('Reading psi(R) at k = {!s}, E = {}\n'.format(K, E_str))
            sys.stdout.flush()

        for line in psi_file:
            line_split = line.split()
            i = int(line_split[0])
            j = int(line_split[1])
            real = float(line_split[2])
            imag = float(line_split[3])
            psi[i, j] = complex(real, imag)
        return psi

    def calculate_current_scan(self, z, V, T, tip_work_func, tip_energy, delta_s=None,
                               fraction=0.025, recalculate=False, vectorised=True,
                               interpolation='cubic', debug=False):
        """Calculate tunnelling current across a plane ie. in constant-height mode.

        Args:
            z (float): z-value of plane in Bohr radii; Uses nearest mesh point to given value
            V (float): Bias voltage
            T (float): Absolute temperature
            tip_work_func (float): Work function of tip
            tip_energy (float): Fermi-level of tip
            delta_s (float, opt): Surface broadening parameter; If None, uses default value
            fraction (float, opt): Fraction of max charge density to use as value for isosurface
            recalculate (bool, opt): Force recalculation, even if already stored.
            vectorised (bool, opt): If true, use NumPy vectorisation.
            interpolation (str, opt): Method of interpolation.
                Possible arguments are 'linear', 'quadratic', 'cubic'.
            debug (bool, opt): Print extra information during runtime.

        Returns:
            array(float): 2D array of current values.
        """

        if delta_s is None:
            delta_s = self.default_delta_s

        if debug:
            sys.stdout.write('Calculating I(R) at V={}V\n'.format(V))
            sys.stdout.flush()

        min_E, max_E = self.bias_to_energy_range(V)

        total_k_weight = 0
        total_energies = 0
        for K in self.bands:
            total_k_weight += K.weight
            for E in self.bands[K]:
                if min_E <= E <= max_E:
                    total_energies += 1
        energies_done = 0

        # Find nearest mesh z-value to given
        z, k = self.get_nearest_mesh_value(z, points=self.real_mesh.shape[2])

        # Initialise meshes
        current = np.zeros(self.real_mesh.shape[:2], dtype=float)
        psi = np.zeros_like(current, dtype=complex)
        elements = current.shape[0]*current.shape[1]

        # Calculate meshes
        ldos = self.get_ldos_grid(min_E, max_E, T, recalculate=recalculate, vectorised=vectorised,
                                  interpolation=interpolation, debug=debug)
        c = self.get_c(ldos, fraction, k, delta_s=delta_s)
        G_conjugate = np.conjugate(self.greens_function_mesh(k, tip_work_func,
                                                             tip_energy, debug=debug))

        # Print debug info
        if debug:
            sys.stdout.write('Calculating grad(G)\n')
            sys.stdout.flush()

        G_conjugate_gradient = self.periodic_gradient(G_conjugate)
        G_centre = G_conjugate.shape[0] / 2

        # Iterate over energies
        for K in self.bands:
            w = K.weight
            for E in self.bands[K]:
                if min_E <= E <= max_E:
                    fd = self.fermi_dirac(E - V, T)
                    if V < 0:
                        fd = 1 - fd

                    # Get propagated wavefunction
                    if not recalculate and os.path.isfile(
                            self.propagated_psi_filename(K, E, T, fraction, z, delta_s=delta_s)):
                        # Read data from file
                        psi = self.read_prop_psi(K, E, T, fraction, z, delta_s=delta_s, debug=debug)
                        read = True
                    else:
                        read = False
                        points_done = 0
                        bars_done = 0

                        # Get unpropagated wavefunction
                        raw_psi = self.get_psi_grid(K, E, recalculate=recalculate,
                                                    vectorised=vectorised,
                                                    interpolation=interpolation, debug=debug)

                        # Print debug info
                        prog = float(energies_done) / total_energies * 100
                        if self.PRINT_RELATIVE_TO_EF:
                            E_str = str(E - self.fermi_level) + ' eV'
                        else:
                            E_str = str(E) + ' eV'
                        debug_str = 'Calculating psi(R) at k = {!s}, E = {}: {:5.1f}%'.format(K,
                                                                                              E_str,
                                                                                              prog)
                        if debug:
                            sys.stdout.write(debug_str)
                            sys.stdout.flush()

                        # Get meshes
                        A = self.get_A_mesh(c, raw_psi)
                        B = self.get_B_mesh(c, raw_psi)

                        # Iterate over tip positions
                        for i, j in np.ndindex(current.shape):
                            # By rolling the tip Green's function, we move the position of the tip
                            G_conjugate_rolled = np.roll(G_conjugate, (i - G_centre), 0)
                            G_conjugate_rolled = np.roll(G_conjugate_rolled, (j - G_centre), 1)[
                                                 :self.real_mesh.shape[0], :self.real_mesh.shape[1]]
                            G_conjugate_gradient_rolled = np.roll(G_conjugate_gradient, (i - G_centre), 0)
                            G_conjugate_gradient_rolled = np.roll(G_conjugate_gradient_rolled, (j - G_centre), 1)[
                                                          :self.real_mesh.shape[0], :self.real_mesh.shape[1]]

                            # Perform volume integral
                            integrand = G_conjugate_rolled * A - self.mesh_dot_product(B, G_conjugate_gradient_rolled)
                            psi[i, j] = np.sum(integrand) * self.grid_spacing ** 3

                            # Update progress bar
                            points_done += 1
                            prog = float(points_done) / elements
                            if debug and prog * self.PROG_BAR_INTERVALS >= bars_done:
                                percent = prog * 100
                                sys.stdout.write('\r')
                                sys.stdout.write(debug_str)
                                sys.stdout.write(' [{:<{}}] {:3.0f}%'.format(self.PROG_BAR_CHARACTER * bars_done,
                                                                             self.PROG_BAR_INTERVALS, percent))
                                sys.stdout.flush()
                                bars_done += 1
                        # Save propagated wavefunction
                        self.write_prop_psi(psi, K, E, T, fraction, z, delta_s=delta_s)

                    # Add to current total
                    current += fd * (w / total_k_weight) * abs(psi)**2

                    energies_done += 1
                    if debug and not read:
                        sys.stdout.write("\n")
                        sys.stdout.flush()
        return current

    def current_filename(self, z, V, T, fraction, delta_s=None):
        """Return standardised filename for relevant current file."""
        
        if delta_s is None:
            delta_s = self.default_delta_s
        return (self.MESH_FOLDER+self.CURRENT_FNAME+self.name+'_'+str(self.grid_spacing)+'_'+str(z)+
                '_'+str(V)+'_'+str(T)+'_'+str(fraction)+'_'+str(delta_s)+self.EXT)

    def write_current(self, current, z, V, T, fraction, delta_s=None, debug=False):
        """Write current to file."""
        
        filename = self.current_filename(z, V, T, fraction, delta_s)

        current_file = safe_open(filename, 'w')
        if debug:
            sys.stdout.write('Writing current grid to {}\n'.format(filename))
        # Iterate over mesh points
        for i, j in np.ndindex(current.shape):
            # If LDOS is non-zero at mesh point, write data to file
            if current[i, j] != 0:
                current_file.write(str(i)+" "+str(j)+" "+str(current[i, j])+'\n')
        current_file.close()

    def read_current(self, z, V, T, fraction, delta_s=None, debug=False):
        """Read current grid from file."""
        
        filename = self.current_filename(z, V, T, fraction, delta_s=delta_s)
        current_file = open(filename, 'r')
        current = np.zeros(self.real_mesh.shape[:2], dtype=float)

        if debug:
            sys.stdout.write('Reading I(R) from {}\n'.format(filename))
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
            sys.stdout.write('I(R) successfully read\n')
            sys.stdout.flush()
        return current

    def get_current_scan(self, z, V, T, tip_work_func, tip_energy, delta_s=None, fraction=0.025,
                         recalculate=False, vectorised=True, interpolation='cubic', debug=False):
        """Get constant-height tunnelling current as 2d array.

        If it is already stored in a file and recalculate is False, then it will read from file.
        Otherwise, it will calculate current.

        Args:
            z (float): z-value of plane in Bohr radii; Uses nearest mesh point to given value.
            V (float): Bias voltage.
            T (float): Absolute temperature.
            tip_work_func (float): Work function of tip.
            tip_energy (float): Fermi-level of tip.
            delta_s (float, opt): Surface broadening parameter; If None, uses default value.
            fraction (float, opt): Fraction of max charge density to use as value for isosurface.
            recalculate (bool, opt): Force recalculation, even if already stored.
            vectorised (bool, opt): If true, use NumPy vectorisation.
            debug (bool, opt): Print extra information during runtime.

        Returns:
            array(float): 2D array of current values.
        """
        
        if delta_s is None:
            delta_s = self.default_delta_s
        if not recalculate and os.path.isfile(self.current_filename(z, V, T, fraction)):
            # Read data from file
            current = self.read_current(z, V, T, fraction, debug=debug)
        else:
            current = self.calculate_current_scan(z, V, T, tip_work_func, tip_energy,
                                                  delta_s=delta_s, fraction=fraction,
                                                  recalculate=recalculate,
                                                  interpolation=interpolation,
                                                  vectorised=vectorised, debug=debug)
            self.write_current(current, z, V, T, fraction)
        return current

    def get_spectrum_th(self, xy, min_V, max_V, sigma, T, fraction, z, delta_s=None, dE=0.005, debug=False):
        """Get spectroscopic data from a list of specific tip positions.

        Args:
            xy (list(list(float)): x-y points of tip in a0; Given as [[x1, y1], [x2, y2], ...].
                Uses nearest mesh point.
            min_V (float): Lower bound for voltage range.
            max_V (float): Upper bound for voltage range.
            sigma (float): State smearing parameter in eV.
            T (float): Absolute temperature in K.
            z (float): z-value of plane in a0; Uses nearest mesh point to given value.
            fraction (float): Fraction of maximum charge density to use as value for isosurface.
            delta_s (float, opt.): Surface broadening parameter; If None, uses default value.
            dE (float, opt.): Energy resolution of data points.
            debug (bool, opt.): Print extra information during runtime.

        Returns:
            array(float): Range of voltage values with resolution dE.
            array(float): Tunnelling conductance values for each tip position.
        """
        
        # Convert voltage range into absolute energies
        min_E = min_V + self.fermi_level
        max_E = max_V + self.fermi_level

        # Get positions of tip
        mesh_positions = []
        mesh_indices = []
        for l in range(len(xy)):
            x, i = self.get_nearest_mesh_value(xy[l][0], points=self.real_mesh.shape[0])
            y, j = self.get_nearest_mesh_value(xy[l][1], points=self.real_mesh.shape[1])
            mesh_positions.append((x, y))
            mesh_indices.append((i, j))

        Es = []
        psis = []
        weights = []
        l = 0
        # Iterate over energies in voltage range
        for K in self.bands:
            for E in self.bands[K]:
                # Assume states 3*sigma away from the edges contribute negligibly
                if min_E - 3*sigma < E < max_E + 3*sigma:
                    Es.append(E)
                    weights.append(K.weight)
                    psi = self.read_prop_psi(K, E, T, fraction, z, delta_s=delta_s, debug=debug)
                    psis.append([])

                    # Get wavefunction for each tip position
                    for m in range(len(mesh_positions)):
                        psis[l].append(abs(psi[mesh_indices[m]])**2)
                    l += 1
        # Convert lists to arrays
        Es = np.array(Es)
        psis = np.array(psis)
        weights = np.array(weights)

        # Generate energy points for plot
        E_range = np.arange(min_E, max_E, dE)
        # Initialise LDOS
        LDOS = np.zeros((E_range.shape[0], len(xy)))

        for u in range(len(E_range)):
            LDOS[u] = np.sum(weights[..., None] * psis *
                             np.exp(-(((E_range[u] - Es[..., None]) / sigma)**2) / 2), axis=0)
        V_range = E_range - self.fermi_level

        return V_range, LDOS

    def get_cits(self, V, T, fraction, sigma, delta_s=None, debug=False):
        """Calculate Current Imaging Tunnelling Spectroscopy scan - UNFINISHED and UNRELIABLE."""
        
        if delta_s is None:
            delta_s = self.default_delta_s

        min_E, max_E = self.bias_to_energy_range(V)

        ldos = self.get_ldos_grid(min_E, max_E, T, debug=debug)
        b = self.broadened_surface(ldos, fraction, self.real_mesh[..., -1, 2], delta_s=delta_s)

        scan = np.zeros(b.shape[:2])

        for K in self.bands:
            w = K.weight
            for E in self.bands[K]:
                if -3*sigma < E - V - self.fermi_level < 3*sigma:
                    exp = np.exp(-(((E - V - self.fermi_level) / sigma)**2) / 2)
                    psi = self.get_psi_grid(K, E, debug=debug)
                    b_psi = b * abs(psi)**2

                    for i, j in np.ndindex(scan.shape):
                        scan[i, j] += w * exp * np.sum(b_psi[i, j])
        scan = scan.reshape(b.shape[:2])

        return scan
