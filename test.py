
from packages import atomic
from packages import io
from packages import cell

import glob

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
Prsr = io.Parser()
Prsr.parseIons(ionFolder, ionFiles)

conqFiles = ['C6H6_SZ']
#print Prsr.getIon('H_SZ_6.5au').sortPAOs()
Prsr.parseConq(conqFolder, conqFiles)

print 'Parsed OK'

atoms = Prsr.atoms

for i in range(1, 2):
	x = atoms[i].x
	y = atoms[i].y
	z = atoms[i].z
	print x, y, z

C = cell.Cell(20.0, 20.0, 20.0, gridSpacing=0.5)

for i in range(1, 13):
	C.addAtom(atoms[i], i)

C.setPsi()
print C.psi

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

