
class Radial:

	def __init__(self, zeta, n, l, r, R, cutoff):
		self.zeta = zeta
		self.n = n
		self.l = l
		self.r = r
		self.R = R
		self.cutoff = cutoff
	
class Ion:

	def __init__(self, name):
		self.name = name
		self.nl = {}
		self.zetas = 1
		self.cutoffs = []
		self.Rads = {}

	def setRadial(self, radial):
		zeta = radial.zeta
		n = radial.n
		l = radial.l
		cutoff = radial.cutoff
		if cutoff not in self.cutoffs:
			self.cutoffs.append(cutoff)
			++self.zetas

		if not self.Rads.has_key(zeta):
			self.Rads[zeta] = {}
		if not self.Rads[zeta].has_key(n):
			self.Rads[zeta][n] = {}
			self.nl[n] = []
		self.Rads[zeta][n][l] = radial
		self.nl[n].append(l)

	#def setRadialFromData(self, n, l, r, R, cutoff):
	#	self.Rads = Radial(n, l, r, R, cutoff)

	def getRadial(self, zeta, n, l):
		return self.Rads[zeta][n][l]

