
from ion import *
import math

class Parser:

	def __init__(self, folder, files):
		self.ionfolder = folder
		self.ionfiles = files
		self.ions = {}

	@property
	def ionfolder(self):
		return self._ionfolder

	def getRadial(self, ionname):
		return self.ions[ionname]

	def getIon(self, ionname):
		return self.ions[ionname]
	
	def parse(self):
		for iname in self.ionfiles:
			# Open file and initialise entry
			f = open(self.ionfolder+iname+'.ion', 'r')

			# Skip preamble and first 9 lines
			line = f.next()
			while not '</preamble>' in line:
				line = f.next()
			for i in range(0, 9):
				line = f.next()

			ion = Ion(iname)
			# Parse PAO data
			line = f.next()
			while line.split()[0] != '#':
				metadata = line.split()
				l = int(metadata[0])
				n = int(metadata[1])
				zeta = int(metadata[2])

				line = f.next()
				metadata = line.split()
				pts = int(metadata[0])
				cutoff = float(metadata[2])

				#self.Rnl[ion][orbitaldata] = [[], [], cutoff]
				r = []
				R = []
				for i in range(0, pts):
					line = f.next()
					x, y = line.split()
					x = float(x)
					y = float(y)
					y = y * math.pow(x, l)
					r.append(x)
					R.append(y)

				R = Radial(zeta, n, l, r, R, cutoff)

					#self.Rnl[ion][orbitaldata][0].append(float(a))
					#self.Rnl[ion][orbitaldata][1].append(float(b)*math.pow(float(a),int(l)))
				ion.setRadial(R)
				line = f.next()
			f.close()
			self.ions[iname] = ion

			