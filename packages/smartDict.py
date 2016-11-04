
class SmartDict(dict):
	"""Extends dict to create nested dicts if invalid key is used rather than throw
	an exception.

	eg. D = SmartDict()
		print D[1][2][4]
		>> {1 : {2 : {4 : {} } } }

		E = SmartDict()
		E[1][3] = "Hello World"
		print E[1][3]
		>> Hello World
	"""

	def __missing__(self, key):
		"""For invalid key, assign empty SmartDict to self[key]"""
		value = self[key] = type(self)()
		return value
