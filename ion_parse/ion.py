
class Radial:

	def __init__(self, n, l, r, R, cutoff):
		self.n = n
		self.l = l
		self.r = r
		self.R = R
		self.cutoff = cutoff

	@property
	def r(self):
		return self.__r

	@property
	def R(self):
		return self.__R
	
	@property
	def n(self):
		return self.__n

	@property
	def l(self):
		return self.__l

	@property
	def cutoff(self):
		return self._cutoff
	
	@r.setter
	def r(self, r):
		self.__r = r

	@R.setter
	def R(self, R):
		self.__R = R

	@n.setter
	def n(self, n):
		self.__n = n

	@l.setter
	def l(self, l):
		self.__l = l

	@cutoff.setter
	def cutoff(self, cutoff):
		self.__cutoff = cutoff

	
class Ion:

	def __init__(self, Name, maxn, maxl):
		self.name = name
		self.maxn = maxn
		self.maxl = maxl

	def setRadial(self, n, l, r, R, cutoff):
		self.Rad = Radial(n, l, r, R, cutoff)









