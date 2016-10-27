
class Radial:

	def __init__(self, zeta, n, l, r, R):
		self.zeta = zeta
		self.n = n
		self.l = l
		self.r = r
		self.R = R
	
class Ion:

	def __init__(self, name):
		self.name = name
		self.nl = {}
		self.zetas = 1
		self.Rads = {}

	def setRadial(self, radial):
		zeta = radial.zeta
		n = radial.n
		l = radial.l

		if zeta > self.zetas:
			self.zetas = zeta


		if not self.Rads.has_key(zeta):
			self.Rads[zeta] = {}
		if not self.Rads[zeta].has_key(n):
			self.Rads[zeta][n] = {}
			self.nl[n] = []
		if not l in self.nl[n]:
			self.nl[n].append(l)

		self.Rads[zeta][n][l] = radial

	#def setRadialFromData(self, n, l, r, R, cutoff):
	#	self.Rads = Radial(n, l, r, R, cutoff)

	def getRadial(self, zeta, n, l):
		return self.Rads[zeta][n][l]

