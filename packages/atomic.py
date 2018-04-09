
from scipy.interpolate import interp1d

from smart_dict import SmartDict
from sph import sph
from vector import Vector

__author__ = 'Johan G. McQuillan'
__email__ = 'johan.mcquillan.13@ucl.ac.uk'


class Radial(object):
    """Evaluates radial part of support function by interpolating data read from .ion file.

    All lengths measured in Bohr radii (a0).
    All energies measured in electron volts (eV).

    Attributes:
        n (int): Principal quantum number; Not used for any calculations
        l (int): Orbital angular momentum quantum number
        zeta (int): Indexes cut-off radii
        radii (list(float)): List of radial distance values as given in .ion file
        radial_function_values (list(float)): List of radial function values for each value of radii
        cutoff (float): Range of radial function, beyond which value of R is 0.0
    """

    interpolators = ['linear', 'quadratic', 'cubic']

    def __init__(self, n, l, zeta, radii, radial_function_values, cutoff):
        """Constructs radial part of a support function.

        Args:
            n (int): Principal quantum number
            l (int): Orbital angular momentum quantum number
            zeta (int): Indexes cut-off radii
            radii (list(float)): List of radial distance values as given in .ion file
            radial_function_values (list(float)): List of radial function values for each value of radii
            cutoff (float): Range of radial function, beyond which value of R is 0.0
        """
        self.zeta = zeta
        self.n = n
        self.l = l
        self.radii = radii
        self.radial_function_values = radial_function_values
        self.cutoff = cutoff

        # Use given data to create interpolation methods
        self.interpolator_linear = interp1d(radii, radial_function_values, kind='linear')
        self.interpolator_quadratic = interp1d(radii, radial_function_values, kind='quadratic')
        self.interpolator_cubic = interp1d(radii, radial_function_values, kind='cubic')

    def get_value(self, distance, interpolation='cubic'):
        """Evaluate radial function value at a given distance using specified interpolation method.

        Args:
            distance (float): Distance to evaluate radial function
            interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'

        Returns:
            float: Radial function value
        """
        if distance > self.cutoff:
            value = 0.0
        else:
            if interpolation == 'linear':
                value = self.interpolator_linear(distance)
            elif interpolation == 'quadratic':
                value = self.interpolator_quadratic(distance)
            elif interpolation == 'cubic':
                value = self.interpolator_cubic(distance)
            else:
                errorString = 'Interpolation method must be either '
                for i in range(len(self.interpolators)):
                    errorString += self.interpolators[i]
                    if i < len(self.interpolators) - 1:
                        errorString += ', '
                raise ValueError(errorString)
        return value


class Ion(object):
    """Represents a certain element represented in a certain basis.

    Stores Radial objects read from .ion file. Used as a template by the Atom sub-class.

    All lengths measured in Bohr radii (a0).
    All energies measured in electron volts (eV).

    Attributes:
        ion_name (string): Name of ion (usually basename of .ion file)
        radials (SmartDict): Radial objects accessed by radials[l][zeta]
    """

    def __init__(self, ion_name, radial_dict=SmartDict()):
        """Construct ion representing an element in a given basis set.
        
        Args:
            ion_name (string): Name of ion (usually from name of .ion file)
            radial_dict (SmartDict, opt.): Radial objects, indexed by radials[zeta][n][l],
                                                Default is empty, Radial objects may be added after instantiation
        """
        self.ion_name = ion_name
        self.radials = radial_dict  # Radial objects; accessed by self.radials[zeta][n][l]

    def sorted_pao(self):
        """Sort pseudo-atomic orbitals into order used by .dat files.

        Returns:
            list: Ordered list of PAO data; Each element is a list containing [l, zeta, m] for the PAO
        """
        pao_list = []

        # Dict keys not necessarily in order
        # Order key list by lowest to highest for each index
        l_list = sorted(self.radials.keys())
        for l in l_list:
            zeta_list = sorted(self.radials[l])
            for zeta in zeta_list:
                for m in range(-l, l + 1):
                    pao_list.append((l, zeta, m))
        return pao_list

    def has_radial(self, l, zeta):
        """Check if Ion has radial object for specified object without creating a dict.

        This encapsulation is required due to autovivification of SmartDict.

        Args:
            l (int): Orbital angular momentum quantum number
            zeta (int): Indexes cut-off radii

        Returns:
            boolean: True if Radial is stored, false if not
        """
        output = False
        if l in self.radials:
            if zeta in self.radials[l]:
                return True
        return output

    def add_radial(self, radial):
        """Adds radial to self.radials.

        Overwrites radial with same metadata (l, zeta).

        Args:
            radial (Radial): Radial object to add
        """
        # Get metadata
        zeta = radial.zeta
        l = radial.l

        # Add Radial
        self.radials[l][zeta] = radial
        self.sorted_pao()

    def get_radial(self, l, zeta):
        """Return Radial object for specified orbital.

        This encapsulation is required due to autovivification of SmartDict.
        Otherwise, if self.radials[l][zeta] is called when there is no entry, it will create an empty entry.

        Args:
            l (int): Orbital angular momentum quantum number
            zeta (int): Indexes cut-off radii

        Returns:
            Radial: Radial object for specified indices
        """
        if self.has_radial(l, zeta):
            return self.radials[l][zeta]
        else:
            raise ValueError("No radial data for " + self.ion_name + ", l = " + str(l) + ", zeta = " + str(zeta))

    def get_radial_value(self, l, zeta, distance, interpolation='cubic'):
        """Use interpolation to evaluate radial function at distance r.

        Encapsulation required due to autovivification of SmartDict.

        Args:
            l (int): Orbital angular momentum quantum number
            zeta (int): Indexes cut-off radii
            distance (float): Radial distance from ion
            interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'

        Returns:
            float: Radial function evaluated at r
        """
        return self.get_radial(l, zeta).get_value(distance, interpolation=interpolation)

    def get_max_cutoff(self):
        """Return the maximum cut-off radius for all stored Radials.

        Beyond the cutoff the radial part of the support function is defined to be 0.
        """
        max_cutoff = 0.0
        for l in self.radials:
            for zeta in self.radials[l]:
                if max_cutoff < self.radials[l][zeta].cutoff:
                    max_cutoff = self.radials[l][zeta].cutoff
        return max_cutoff


class Atom(Ion):
    """Stores information on an atom, extending Ion to include atomic position and support function coefficients.

    All lengths measured in Bohr radii (a0).
    All energies measured in electron volts (eV).

    Attributes:
        ion_name (string): Name of ion (usually from name of .ion file)
        radials (SmartDict): Radial objects accessed by radials[l][zeta], where all indices are int
        atom_pos (Vector): Vector for atom position relative to simulation cell
        bands (SmartDict): SmartDict of complex basis coefficients;
                            Accessed by bands[K][E][l][zeta][m], where K is k-point Vector object, and E is energy
    """

    def __init__(self, ion_name, atom_position, radials=SmartDict()):
        """Constructor for atom.

        Args:
            ion_name (string): Name of ion (usually from name of .ion file)
            atom_position (Vector): 3D Cartesian real space vector for atom position
            radials (SmartDict, opt.): Radial objects, indexed by radials[zeta][n][l],
                                            Default is empty, Radial objects may be added after instantiation
        """
        Ion.__init__(self, ion_name, radials)
        self.atom_pos = atom_position
        self.bands = SmartDict()

    def set_ion(self, ion):
        """Copy all attributes from an Ion to this Atom."""
        self.ion_name = ion.ion_name
        self.radials = ion.radials
        self.sorted_pao()

    def within_cutoff_relative(self, relative_position, l=None, zeta=None):
        """Return true if within cut-off region of specified orbital, or if none specified, maximum cut-off region

        Args:
            relative_position (Vector): Vector relative to atom position
            l (int, opt.): Orbital angular momentum quantum number; to check specific radial, needs zeta
            zeta (int, opt.): Zeta index; to check specific radial, needs l

        Returns:
            bool: True if within cutoff radius
        """
        output = False
        distance = abs(relative_position)
        if l is None or zeta is None:
            # No orbital specified
            # Check for max cut-off
            if distance <= self.get_max_cutoff():
                output = True
        else:
            if self.has_radial(l, zeta):
                if distance <= self.get_radial(l, zeta).cutoff:
                    output = True
        return output

    def within_cutoff(self, position, l=None, zeta=None):
        """Return true if within cutoff region.

        Args:
            position (Vector): Vector relative to simulation cell origin
            l (int, opt.): Orbital angular momentum quantum number; to check specific radial, needs zeta
            zeta (int, opt.): Indexes cut-off radii

        Returns:
            bool: True if within cut-off region
        """
        relative_position = position - self.atom_pos
        return self.within_cutoff_relative(relative_position, l, zeta)

    def has_coefficient(self, K, E, l, zeta, m):
        """Check if atom stores coefficient for given orbital.

        Encapsulation required due to autovivification of SmartDict.

        Args:
            K (Vector): K-space Vector
            E (float): Band energy
            l (int): Orbital angular momentum quantum number
            zeta (int): Indexes cut-off radii
            m (int): Magnetic quantum number

        Returns:
            boolean: True if coefficient is stored, false if not
        """
        output = False
        if K in self.bands:
            if E in self.bands[K]:
                if l in self.bands[K][E]:
                    if zeta in self.bands[K][E][l]:
                        if m in self.bands[K][E][l][zeta]:
                            output = True
        return output

    def add_coefficient(self, K, E, PAO, coefficient):
        """Add a complex coefficient to self.bands

        Args:
            K (Vector): K-point Vector
            E (float): Band energy
            PAO (int): Index of PAO as given in .dat file
            coefficient (complex): Coefficient of PAO
        """
        PAOdata = self.sorted_pao()[PAO - 1]
        l = PAOdata[0]
        zeta = PAOdata[1]
        m = PAOdata[2]

        self.bands[K][E][l][zeta][m] = coefficient

    def get_coefficient(self, K, E, l, zeta, m):
        """Return complex coefficient for given orbital.

        Args:
            K (Vector): K-point Vector
            E (float): Band energy
            l (int): Orbital angular momentum quantum number
            zeta (int): Indexes cut-off radii
            m (int): Magnetic quantum number

        Returns:
            complex: Coefficient for given orbital
        """
        output = 0.0
        if self.has_coefficient(K, E, l, zeta, m):
            output = self.bands[K][E][l][zeta][m]
        return output

    def get_radial_value_relative(self, l, zeta, relative_position, interpolation='cubic'):
        """Evaluate radial part of support function.

        Args:
            l (int): Orbital angular momentum quantum number
            zeta (int): Indexes cut-off radii
            relative_position (Vector): Vector relative to atom
            interpolation (string, opt.): Method of interpolation; possible arguments are 'linear', 'quadratic', 'cubic'

        Returns:
            float: Value of radial part of support function
        """
        R = 0.0
        if self.has_radial(l, zeta):
            distance = abs(relative_position)
            R = self.get_radial_value(l, zeta, distance, interpolation=interpolation)
        return R

    def get_sph(self, l, m, position):
        """Evaluate real spherical harmonic with atom as origin.

        Args:
            l (int): Orbital angular momentum quantum number
            m (int): Magnetic quantum number
            position (Vector): Vector relative to simulation cell origin

        Returns:
            float: Value of real spherical harmonic
        """
        relative_position = position - self.atom_pos
        return sph(l, m, relative_position)
