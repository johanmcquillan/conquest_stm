import math

from .. import atomic
from ..smartDict import SmartDict

class Parser(object):

	"""Parses and stores data from input files.
	Currently only works for .ion files."""

	def __init__(self, ionFolder, ionFiles, conqFolder, conqFiles):
		"""Constructor for Parser.

		Args:
			ionFolder (string): Relative folder path to .ion files
			ionFiles ([string]): Ion filenames, excluding '.ion'
			conqFolder (string): Relative folder path to .dat and Conquest_out files
			conqFiles ([string]): Conquest ouptut filenames, excluding '.dat'
		"""

		self.ions = {}
		self.atoms = SmartDict()
		self.fermiLevels = {}

		self.ionFiles = ionFiles
		self.ionFolder = ionFolder
		self.conqFiles = conqFiles
		self.conqFolder = conqFolder

	def parseIons(self):
		"""Parse data from ionFiles to Ion objects and store in self.ions indexed by ionFile name."""

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

	def parseConq(self, tolerance=0.0, debug=False):
		"""Parse data from conqFiles to Atom objects and store in self.ions indexed by conqFile name."""

		# Open Conquest_out file
		for conq in self.conqFiles:
			Fconq = open(self.conqFolder+conq)

			# Skip opening lines
			line = Fconq.next()
			while len(line.split()) != 8:
				line = Fconq.next()

			# Read atomic positions
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

				for atomKey in atomData:
					atomDataList = atomData[atomKey]
					if atomDataList[3] == ionType:
						x, y, z = atomDataList[:3]
						self.atoms[conq][atomKey] = atomic.Atom(ionName, x, y, z)
						self.atoms[conq][atomKey].setIon(self.ions[ionName])
				Fconq.next()
				line = Fconq.next()

			Fconq.close()

			# Open corresponding .dat file
			Fcoeff = open(self.conqFolder+conq+'.dat')
			line = Fcoeff.next()

			endOfFile = False
			while not endOfFile:
				if '#Kpoint' in line:
					line = Fcoeff.next()
					data = line.split()
					k1 = float(data[0])
					k2 = float(data[1])
					k3 = float(data[2])
					kpoint = [k1, k2, k3]
					try:
						line = Fcoeff.next()
						while '#Kpoint' not in line:
							data = line.split()
							bandN = int(data[0])
							bandE = float(data[1])

							# If a tolerance is set, add coefficients to existing band
							if tolerance != 0.0:
								foundExistingBand = False
								addedAtomKeys = self.atoms[conq].keys()
								i = 0
								while not foundExistingBand and i < len(addedAtomKeys):
									addedAtomKey = addedAtomKeys[i]
									addedAtom = self.atoms[conq][addedAtomKey]
									existingBands = addedAtom.bands.keys()
									j = 0
									while not foundExistingBand and j < len(existingBands):
										E = existingBands[j]
										if abs(E - bandE) < tolerance:
											if debug:
												print 'Combining band '+str(bandE)+' with '+str(E)
											bandE = E
											foundExistingBand = True

										j = j + 1
									i = i + 1

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
								self.atoms[conq][a].addCoeff(bandE, PAO, complexCoeff, combine=True)
								line = Fcoeff.next()

					except StopIteration:
						endOfFile = True
				else:
					# Check if end of file
					try:
						line = Fcoeff.next()
					except StopIteration:
						endOfFile = True
			Fcoeff.close()

			FDoS = open(self.conqFolder+conq+'.dos')
			line = FDoS.next()
			FDoS.close()

			data = line.split()
			fermiLevel = float(data[2])

			self.fermiLevels[conq] = fermiLevel
