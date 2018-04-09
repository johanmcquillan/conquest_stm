
import numpy as np

__author__ = 'Johan G. McQuillan'
__email__ = 'johan.mcquillan.13@ucl.ac.uk'

# Normalisation factors for spherical harmonics
# Numbers in variable names refer to l and absolute value of m
SPH0 = 1.0/2.0 * np.sqrt(1.0 / np.pi)

SPH1 = 1.0/2.0 * np.sqrt(3.0 / np.pi)

SPH2 = 1.0/2.0 * np.sqrt(15.0 / np.pi)

SPH33 = 1.0/4.0 * np.sqrt(35.0 / (2.0 * np.pi))
SPH32 = 1.0/2.0 * np.sqrt(105.0 / np.pi)
SPH31 = 1.0/4.0 * np.sqrt(21.0 / (2.0 * np.pi))
SPH30 = 1.0/4.0 * np.sqrt(7.0 / np.pi)

SPH44 = 3.0/4.0 * np.sqrt(35.0 / np.pi)
SPH43 = 3.0/4.0 * np.sqrt(35.0 / (2*np.pi))
SPH42 = 3.0/4.0 * np.sqrt(5.0 / np.pi)
SPH41 = 3.0/4.0 * np.sqrt(5.0 / (2*np.pi))
SPH40 = 3.0/16.0 * np.sqrt(1.0 / np.pi)


def sph(l, m, vector):
    """Return value of real spherical harmonic using Cartesian coordinates.

    Valid indices are l > 0; -l < m < +l.
    Currently, only l <= 4 is supported.

    Args:
        l (int): Degree; Orbital angular momentum quantum number; Must be => 0
        m (int): Order; Azimuthal quantum number; must be within -l >= m >= l
        vector (Vector): 3D real space Cartesian vector
    Returns:
        harmonic (float): Spherical harmonic for given l and m, evaluated at (x, y, z)
    """

    harmonic = 0.0

    # Normalise coordinates to unit magnitude
    normalised_vector = vector.normalise()
    x = normalised_vector.x
    y = normalised_vector.y
    z = normalised_vector.z

    # Choose correct form of spherical harmonic depending on l and m
    if l == 0:
        harmonic = SPH0

    elif l == 1:
        if m == -1:
            harmonic = SPH1 * y
        elif m == 0:
            harmonic = SPH1 * z
        elif m == 1:
            harmonic = SPH1 * x

    elif l == 2:
        if m == -2:
            harmonic = SPH2 * x*y
        elif m == -1:
            harmonic = SPH2 * y*z
        elif m == 0:
            harmonic = SPH2 / (2.0 * np.sqrt(3.0)) * (3*z**2 - 1)
        elif m == 1:
            harmonic = SPH2 * x*z
        elif m == 2:
            harmonic = SPH2 / 2.0 * (x**2 - y**2)

    elif l == 3:
        if m == -3:
            harmonic = SPH33 * y*(3*x**2 - y**2)
        elif m == -2:
            harmonic = SPH32 * x*y*z
        elif m == -1:
            harmonic = SPH31 * y*(5*z**2 - 1)
        elif m == 0:
            harmonic = SPH30 * z*(5*z**2 - 3)
        elif m == 1:
            harmonic = SPH31 * x*(5*z**2 - 1)
        elif m == 2:
            harmonic = SPH32 * z*(x**2 - y**2)
        elif m == 3:
            harmonic = SPH33 * x*(x**2 - 3*y**2)

    elif l == 4:
        if m == -4:
            harmonic = SPH44 * x*y*(x**2 - y**2)
        elif m == -3:
            harmonic = SPH43 * y*z*(3*x**2 - y**2)
        elif m == -2:
            harmonic = SPH42 * x*y*(7*z**2 - 1)
        elif m == -1:
            harmonic = SPH41 * y*z*(7*z**2 - 3)
        elif m == 0:
            harmonic = SPH40 * (35*z**4 - 30*z**2 + 3)
        elif m == 1:
            harmonic = SPH41 * x*z*(7*z**2 - 3)
        elif m == 2:
            harmonic = SPH42/2.0 * (x**2 - y**2)*(7*z**2 - 1)
        elif m == 3:
            harmonic = SPH43 * x*z*(x**2 - 3*y**2)
        elif m == 4:
            harmonic = SPH44/4.0 * (x**2*(x**2 - 3*y**2) - y**2*(3*x**2 - y**2))

    return harmonic
