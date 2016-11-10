import glob

from packages import atomic
from packages import io
from packages import cell

ionFolder = 'ions/'
ionFilesRaw = glob.glob(ionFolder+'*.ion')
ionFiles = []
tempName = ''
for ionNameRaw in ionFilesRaw:
	tempName = ionNameRaw[len(ionFolder):len(ionNameRaw)-4]
	ionFiles.append(tempName)

conqFolder = 'conquest/'
conqFilesRaw = glob.glob(conqFolder+'*.dat')
conqFiles = []
for atomNameRaw in conqFilesRaw:
	tempName = atomNameRaw[len(conqFolder):len(atomNameRaw)-4]
	conqFiles.append(tempName)

#ionFiles = ['C_TZTP_6.5_4.5_2.5au']

# Initialise parser and get data
Prsr = io.Parser(ionFolder, ionFiles, conqFolder, ['C6H6_SZ', 'C6H6_DZDP', 'C6H6_TZTP'])
Prsr.parseIons()
Prsr.parseConq(tolerance=0.0, debug=True)

print 'Parsed OK'

atoms = Prsr.atoms
fl = Prsr.fermiLevels
molecule = 'C6H6_TZTP'

C = cell.Cell('C6H6', fl[molecule], 25.0, 20.0, 25.0, gridSpacing=0.2)

for a in atoms[molecule]:
	C.addAtom(atoms[molecule][a], a)

Pltr = io.Plotter('Test')

for band in range(0, len(C.bands)):
	E = C.bands[band]
	print band, E
	#C.setPsi(E=E, debug=True)
	for axis in ['z']:
		if E > -1.0 and E < 1.0:
			Pltr.plotPsiCrossSec(molecule, C, band, axis, 0, 25, label=str(E), tolerance=0.0, printStatus=True, debug=True)




# Put Ion objects into a dict indexed by name in ionfiles
# ions = {}
# for ion in ionFiles:
# 	ions[ion] = Prsr.getIon(ion)

# Pltr = IO.Plotter('Rnl', ions)
# Pltr.plotRadials()

# I = ions['C_TZTP_6.5_4.5_2.5au']
# #I.plotSPH(2, 1)
# for n in I.nl.keys():
# 	for l in I.nl[n]:
# 		for m in range(-l, l+1):
# 			for e in ['x','y','z']:
# 				I.plotBasis(1, n, l, m, e, step=0.1)
# 				#continue

# for l in range(0, 4):
# 	for m in range(-l, l+1):
# 		for e in ['x', 'y', 'z']:
# 			#Ion.plotSPH(l, m, e)
# 			continue

# # # Plot radial functions to 'Rnl_radials.pdf'
# # Pltr = Plotter('Rnl', ions)
# # Pltr.plotRadials()
