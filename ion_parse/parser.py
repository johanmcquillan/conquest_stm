from ion import *
import math

class Parser:
	"""Parses and stores data from input files.
	Currently only works for .ion files."""

	def __init__(self, ionfolder, ionfiles):
		self.ionfolder = ionfolder
		self.ionfiles = ionfiles
		self.ions = {}

	def getIon(self, iname):
		return self.ions[iname]

	def parseIons(self):
		"""Parse data from ion files to Ion objects and 
		store in self.ions indexed by ionfile name."""
		for iname in self.ionfiles:
			# Open .ion and initialise entry
			f = open(self.ionfolder+iname+'.ion', 'r')

			# Skip preamble and first 9 lines
			line = f.next()
			while not '</preamble>' in line:
				line = f.next()
			for i in range(0, 9):
				line = f.next()

			# Create empty Ion and fill with radial data
			ion = Ion(iname)
			line = f.next()
			while line.split()[0] != '#':
				# Read quantum numbers and zeta index
				metadata = line.split()
				l = int(metadata[0])
				n = int(metadata[1])
				zeta = int(metadata[2])

				# Get number of points for radial and cutoff radius
				line = f.next()
				metadata = line.split()
				pts = int(metadata[0])
				cutoff = float(metadata[2])

				# Initialise R(r) function data
				r = []
				R = []
				# Read data into Radial object and add to Ion
				for i in range(0, pts):
					line = f.next()
					x, y = line.split()
					x = float(x)
					y = float(y)
					y = y * math.pow(x, l)
					r.append(x)
					R.append(y)

				R = Radial(zeta, n, l, r, R)
				ion.addRadial(R)
				line = f.next()

			f.close()
			self.ions[iname] = ion

