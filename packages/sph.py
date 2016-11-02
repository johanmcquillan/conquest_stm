import numpy as np
import math

# Normalisation factors for spherical harmonics
# Numbers in variable names refer to l and absolute value of m
sph00 = 1.0/2.0 * np.sqrt(1.0   / math.pi)

sph11 = 1.0/2.0 * np.sqrt(3.0   / (2.0*math.pi))
sph10 = 1.0/2.0 * np.sqrt(3.0   / math.pi)

sph22 = 1.0/4.0 * np.sqrt(15.0  / (2.0*math.pi)) 
sph21 = 1.0/2.0 * np.sqrt(15.0  / (2.0*math.pi))
sph20 = 1.0/4.0 * np.sqrt(5.0   / math.pi)

sph33 = 1.0/8.0 * np.sqrt(35.0  / math.pi)
sph32 = 1.0/4.0 * np.sqrt(105.0 / (2.0*math.pi))
sph31 = 1.0/8.0 * np.sqrt(21.0  / math.pi)
sph30 = 1.0/4.0 * np.sqrt(7.0   / math.pi)

def sph(l, m, x, y, z):
		"""Return cartesian spherical harmonic.
		Valid indices are l > 0; -l < m < +l.
		Currently, only l <= 3 is supported.
		
		Input:
		l:			Orbital angular momentum quantum number; Must be => 0
		m:			Azimuthal quantum number; must be within -l >= m >= l
		x, y, z:	Cartesian coordinates

		Output:
		harm:		Spherical harmonic for given l and m, evaluated at (x, y, z)"""

		# Normalise coordinates to unit magnitude
		magnitude = np.sqrt(x**2 + y**2 + z**2)
		if magnitude != 0:
			x = x / magnitude
			y = y / magnitude
			z = z / magnitude

		# Choose correct form of spherical harmonic depending on l and m
		harm = 0.0
		if l == 0:
			harm =      sph00
		# l = 0 is the only harmonic with non-zero value at origin
		elif magnitude == 0:
			harm = 		0.0
		elif l == 1:
			if m == 1:
				harm =  sph11 * x
			elif m == -1:
				harm =  sph11 * y
			elif m == 0:
				harm =  sph10 * z
		elif l == 2:
			if m == 2:
				harm =  sph22 * (x**2 - y**2)
			elif m == -2:
				harm =  sph22 * (x**2 - y**2)
			elif m == 1:
				harm =  sph21 * x*z
			elif m == -1:
				harm = -sph21 * x*z
			elif m == 0:
				harm =  sph20 * (3*z**2 - 1)
		elif l == 3:
			if m == 3:
				harm =  -sph33 * (x**3 - 3*x*y**2)
			elif m == -3:
				harm =   sph33 * (x**3 - 3*x*y**2)
			elif m == 2:
				harm = 	 sph32 * (x**2 - y**2)*z
			elif m == -2:
				harm =   sph32 * (x**2 - y**2)*z
			elif m == 1:
				harm =  -sph31 * (5*z**2 - 1)*x
			elif m == -1:
				harm =   sph31 * (5*z**2 - 1)*x
			elif m == 0:
				harm =   sph30 * (5*z**2 - 3)*z
		return harm