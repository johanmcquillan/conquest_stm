
class SmartDict(dict):
	"""Implements Perl's autovivification for dicts.

	Extends dict to create nested dicts if invalid key is used rather than throw
	an exception.

	eg. D = SmartDict()
		print D[1][2][4]
		>> {1 : {2 : {4 : {} } } }

		E = SmartDict()
		E[1][3] = "Hello World"
		print E[1][3]
		>> Hello World
	"""

	def __init__(self):
		dict.__init__(self)
		self.locked = False

	def __missing__(self, key):
		"""If not locked and invalid key, assign empty SmartDict to self[key]"""
		if self.locked:
			raise ValueError(key)
		else:
			value = self[key] = type(self)()
			return value

	def lock(self):
		"""Prevent autovivification"""
		self.locked = True
		for key in self:
			if type(self[key]) is type(self):
				self[key].lock()

	def unlock(self):
		"""Allow autovivification"""
		self.locked = False
		for key in self:
			if type(self[key]) is type(self):
				self[key].unlock()
