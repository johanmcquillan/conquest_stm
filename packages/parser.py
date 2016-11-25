import math

import atomic
from cell import Cell
from smartDict import SmartDict

HA_TO_EV = 0.03674932  # Factor to convert Hartrees to electron volts


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
		self.electrons = {}
		self.cellDimensions = {}

		self.ionFiles = ionFiles
		self.ionFolder = ionFolder
		self.conqFiles = conqFiles
		self.conqFolder = conqFolder

	def parseIons(self):
		"""Parse data from ionFiles to Ion objects and store in self.ions indexed by ionFile name."""

		for ionName in self.ionFiles:
			# Open .ion and initialise entry
			Fion = open(self.ionFolder+ionName+'.ion', 'r')

			# Skip preamble and first 9 lines
			line = Fion.next()
			while '</preamble>' not in line:
				line = Fion.next()
			for i in range(0, 9):
				Fion.next()

			# Create empty dict of radial objects
			radialDict = SmartDict()
			line = Fion.next()
			while line.split()[0] != '#':
				# Read quantum numbers and zeta index
				metadata = line.split()
				l = int(metadata[0])
				n = int(metadata[1])
				zeta = int(metadata[2])

				# Get number of points for radial and cutoff radius
				line = Fion.next()
				metadata = line.split()
				pts = int(metadata[0])
				cutoff = float(metadata[2])

				# Initialise R(r) function data
				r = []
				R = []
				# Read data into Radial object and add to Ion
				for i in range(0, pts):
					line = Fion.next()
					x, y = line.split()
					x = float(x)
					y = float(y) * math.pow(x, l)
					r.append(x)
					R.append(y)

				# Create Radial object and store in dict
				radialDict[l][zeta] = atomic.Radial(n, l, zeta, r, R, cutoff)
				line = Fion.next()
			Fion.close()

			# Create Ion with radials and add to dict
			self.ions[ionName] = atomic.Ion(ionName, radialDict)
			del radialDict

	def parseConquestOutput(self, tolerance=0.0, debug=False):
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

			while "The simulation box has the following dimensions" not in line:
				line = Fconq.next()
			line = Fconq.next()
			rawData = line.split()
			cellLengthX = float(rawData[2])
			cellLengthY = float(rawData[5])
			cellLengthZ = float(rawData[8])
			self.cellDimensions[conq] = [cellLengthX, cellLengthY, cellLengthZ]

			# 
			while '------------------------------------------------------------------' not in line:
				line = Fconq.next()
			Fconq.next()
			Fconq.next()
			line = Fconq.next()

			# Get ion species
			while '------------------------------------------------------------------' not in line:
				rawData = line.split()

				ionType = int(rawData[0])
				ionName = rawData[1]
				numberOfElectrons = int(rawData[5])

				for atomKey in atomData:
					atomDataList = atomData[atomKey]
					if atomDataList[3] == ionType:
						x, y, z = atomDataList[:3]
						self.atoms[conq][atomKey] = atomic.Atom(ionName, x, y, z)
						self.atoms[conq][atomKey].setIon(self.ions[ionName])
						if conq not in self.electrons:
							self.electrons[conq] = 0
						self.electrons[conq] += numberOfElectrons
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
					Kx = float(data[0])
					Ky = float(data[1])
					Kz = float(data[2])
					try:
						line = Fcoeff.next()
						while '#Kpoint' not in line:
							data = line.split()
							bandN = int(data[0])
							bandE = float(data[1])*HA_TO_EV

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

										j += 1
									i += 1

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
								self.atoms[conq][a].addCoeff(Kx, Ky, Kz, bandE, PAO, complexCoeff)
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
			fermiLevel = float(data[2])*HA_TO_EV

			self.fermiLevels[conq] = fermiLevel

	def getCell(self, conq, gridSpacing=0.5):

		electrons = self.electrons[conq]
		fermi = self.fermiLevels[conq]
		x = self.cellDimensions[conq][0]
		y = self.cellDimensions[conq][1]
		z = self.cellDimensions[conq][2]
		
		C = Cell(conq, fermi, electrons, x, y, z, gridSpacing=gridSpacing)

		for atomKey in self.atoms[conq]:
			C.addAtom(self.atoms[conq][atomKey], atomKey)

		return C

