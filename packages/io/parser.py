
from .. import atomic as atm
import math

class Parser:

	"""Parses and stores data from input files.
	Currently only works for .ion files."""

	def __init__(self):
		self.ions = {}
		self.atoms = {}

	def getIon(self, ionName):
		return self.ions[ionName]

	def parseIons(self, ionFolder, ionFiles):
		"""Parse data from ion files to Ion objects and 
		store in self.ions indexed by ionFile name."""
		self.ionFolder = ionFolder
		self.ionFiles = ionFiles

		for ionName in self.ionFiles:
			# Open .ion and initialise entry
			f = open(self.ionFolder+ionName+'.ion', 'r')

			# Skip preamble and first 9 lines
			line = f.next()
			while not '</preamble>' in line:
				line = f.next()
			for i in range(0, 9):
				line = f.next()

			# Create empty Ion and fill with radial data
			ion = atm.Ion(ionName)
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

				R = atm.Radial(zeta, n, l, r, R, cutoff)
				ion.addRadial(R)
				line = f.next()

			f.close()
			self.ions[ionName] = ion

	# def parseAtoms(self, atomFolder, atomFiles):
	# 	"""Parse data from .dat files to Atom objects and 
	# 	store in self.atoms indexed by atomFile name.
	# 	UNFINISHED"""
	# 	self.atomFolder = atomFolder
	# 	self.atomFiles = atomFiles

	# 	for atomName in self.atomFiles:
	# 		f = open(self.atomFolder+atomName+'.dat')

	# 		line = f.next()
	# 		line = f.next()
	# 		data = line.split()
	# 		kpoint = [float(data[0]), float(data[1]), float(data[2])]

	# 		line = f.next()
	# 		data = line.split()
	# 		band = [int(data[0]), float(data[1])]

	# 		line = f.next()
	# 		data = line.split()
	# 		a = int(data[0])
	# 		PAO = int(data[1])
	# 		coeffString = data[2]
	# 		coeffString = coeffString.replace('(', '')
	# 		coeffString = coeffString.replace(')', '')
	# 		complexString = coeffString.split(',')
	# 		complexCoeff = [float(complexString[0]), float(complexString[1])]

	# 		print complexCoeff
	# 		f.close()

	def parseConq(self, conqFolder, conqFiles):
		"""Parse data about atoms store in self.atoms 
		indexes by atomFile name."""
		self.conqFolder = conqFolder
		self.conqFiles = conqFiles

		for conq in self.conqFiles:
			Fconq = open(self.conqFolder+conq)

			line = Fconq.next()
			while len(line.split()) != 8:
				line = Fconq.next()

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

			while line.split() != ['------------------------------------------------------------------']:
				line = Fconq.next()

			Fconq.next()
			Fconq.next()
			line = Fconq.next()

			while line.split() != ['------------------------------------------------------------------']:
				rawData = line.split()

				a = int(rawData[0])
				ionName = rawData[1]

				Fconq.next()
				line = Fconq.next()

			Fconq.close()

			Fcoeff = open(self.conqFolder+conq+'.dat')

			line = Fcoeff.next()
			line = Fcoeff.next()
			data = line.split()
			kpoint = [float(data[0]), float(data[1]), float(data[2])]

			line = Fcoeff.next()
			data = line.split()
			band = [int(data[0]), float(data[1])]

			line = Fcoeff.next()
			data = line.split()
			a = int(data[0])
			PAO = int(data[1])
			coeffString = data[2]
			coeffString = coeffString.replace('(', '')
			coeffString = coeffString.replace(')', '')
			complexString = coeffString.split(',')
			complexCoeff = [float(complexString[0]), float(complexString[1])]

			print complexCoeff
			Fcoeff.close()











