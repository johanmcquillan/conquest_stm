import math

from .. import atomic
from ..smartDict import SmartDict

class Parser(object):

	"""Parses and stores data from input files.
	Currently only works for .ion files."""

	def __init__(self, ionFolder, ionFiles, conqFolder, conqFiles):
		self.ions = {}
		self.atoms = {}

		self.ionFiles = ionFiles
		self.ionFolder = ionFolder
		self.conqFiles = conqFiles
		self.conqFolder = conqFolder

	def parseIons(self):
		"""Parse data from ion files to Ion objects and
		store in self.ions indexed by ionFile name."""

		for ionName in self.ionFiles:
			# Open .ion and initialise entry
			f = open(self.ionFolder+ionName+'.ion', 'r')

			# Skip preamble and first 9 lines
			line = f.next()
			while not '</preamble>' in line:
				line = f.next()
			for i in range(0, 9):
				line = f.next()

			# Create empty dict of radial objects
			radialDict = SmartDict()
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

				# Create Radial object and store in dict
				radialDict[zeta][n][l] = atomic.Radial(zeta, n, l, r, R, cutoff)
				line = f.next()
			f.close()

			# Create Ion with radials and add to dict
			self.ions[ionName] = atomic.Ion(ionName, radialDict)
			del radialDict

	def parseConq(self):
		"""Parse data about atoms store in self.atoms
		indexes by atomFile name."""

		# Open Conquest_out file
		for conq in self.conqFiles:
			Fconq = open(self.conqFolder+conq)

			# Skip opening lines
			line = Fconq.next()
			while len(line.split()) != 8:
				line = Fconq.next()

			# Get atomic positions
			atomData = {}
			while len(line.split()) == 8 and line.split()[0].isdigit():
				rawData = line.split()

				atomIndex = int(rawData[0])
				x = float(rawData[1])
				y = float(rawData[2])
				z = float(rawData[3])
				ionType = int(rawData[4])

				atomData[atomIndex] = [x, y, z, ionType]

				line = Fconq.next()

			# Skip lines until more atom data
			while line.split() != ['------------------------------------------------------------------']:
				line = Fconq.next()
			Fconq.next()
			Fconq.next()
			line = Fconq.next()

			# Get ion species
			while line.split() != ['------------------------------------------------------------------']:
				rawData = line.split()

				ionType = int(rawData[0])
				ionName = rawData[1]

				for atomKey, atomDataList in atomData:
					if atomDataList[3] == ionType:
						x, y, z = atomDataList[:3]
						self.atoms[atomKey] = atomic.Atom(ionName, x, y, z)
						self.atoms[atomKey].setIon(self.ions[ionName])
				Fconq.next()
				line = Fconq.next()

			Fconq.close()

			# Open corresponding .dat file
			Fcoeff = open(self.conqFolder+conq+'.dat')
			line = Fcoeff.next()

			gotCoeffs = False
			while not gotCoeffs:
				line = Fcoeff.next()
				if '#Kpoint' not in line and len(line.split()) != 3:
					data = line.split()
					bandN = int(data[0])
					bandE = float(data[1])

					line = Fcoeff.next()
					while len(line.split()) > 2:
						data = line.split()
						a = int(data[0])
						PAO = int(data[1])
						coeffString = data[2]
						coeffString = coeffString.replace('(', '')
						coeffString = coeffString.replace(')', '')
						complexString = coeffString.split(',')
						complexCoeff = complex(float(complexString[0]), float(complexString[1]))
						self.atoms[a].addCoeff(PAO, complexCoeff)
						line = Fcoeff.next()
					gotCoeffs = True
				else:
					line = Fcoeff.next()
			Fcoeff.close()
