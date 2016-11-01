import numpy as np
import math

def sph(cls, l, m, x, y, z):
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
		# If at origin, return 0

		magnitude = np.sqrt(x**2 + y**2 + z**2)
		if magnitude != 0:
			x = x / magnitude
			y = y / magnitude
			z = z / magnitude
		else:
			return 0.0

		# Choose correct form of spherical harmonic depending on l and m
		harm = 0.0
		if l == 0:
			harm =      cls.sph00
		elif l == 1:
			if m == 1:
				harm =  cls.sph11 * x
			elif m == -1:
				harm =  cls.sph11 * y
			elif m == 0:
				harm =  cls.sph10 * z
		elif l == 2:
			if m == 2:
				harm =  cls.sph22 * (x**2 - y**2)
			elif m == -2:
				harm =  cls.sph22 * (x**2 - y**2)
			elif m == 1:
				harm =  cls.sph21 * x*z
			elif m == -1:
				harm = -cls.sph21 * x*z
			elif m == 0:
				harm =  cls.sph20 * (3*z**2 - 1)
		elif l == 3:
			if m == 3:
				harm =  -cls.sph33 * (x**3 - 3*x*y**2)
			elif m == -3:
				harm =   cls.sph33 * (x**3 - 3*x*y**2)
			elif m == 2:
				harm = 	 cls.sph32 * (x**2 - y**2)*z
			elif m == -2:
				harm =   cls.sph32 * (x**2 - y**2)*z
			elif m == 1:
				harm =  -cls.sph31 * (5*z**2 - 1)*x
			elif m == -1:
				harm =   cls.sph31 * (5*z**2 - 1)*x
			elif m == 0:
				harm =   cls.sph30 * (5*z**2 - 3)*z
		return harm