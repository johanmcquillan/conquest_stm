
class Radial:
	"""Stores the radial part of basis function and metadata, ie.
	quantum numbers (n and l) and zeta index"""

	def __init__(self, zeta, n, l, r, R):
		self.zeta = zeta
		self.n = n
		self.l = l
		self.r = r
		self.R = R
	
class Ion:
	"""Stores Radial objects"""

	def __init__(self, name):
		self.name = name
		self.nl = {}
		self.zetas = 1 # 1 = SZ; 2 = DZ; 3 = TZ
		self.Rads = {} # Radial objects; accessed by self.Rads[zeta][n][l]

	def addRadial(self, radial):
		"""Adds Radial to self.Rads. Overwrites radial with same metadata"""
		# Get metadata
		zeta = radial.zeta
		n = radial.n
		l = radial.l

		# Set max value of zeta
		if zeta > self.zetas:
			self.zetas = zeta

		# Initialise dict entry
		if not self.Rads.has_key(zeta):
			self.Rads[zeta] = {}
		if not self.Rads[zeta].has_key(n):
			self.Rads[zeta][n] = {}
			self.nl[n] = []
		if not l in self.nl[n]:
			self.nl[n].append(l)

		# Add Radial
		self.Rads[zeta][n][l] = radial

	def getRadial(self, zeta, n, l):
		return self.Rads[zeta][n][l]

