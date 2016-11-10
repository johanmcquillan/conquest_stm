
from packages import io, cell
import matplotlib.pyplot as plt

prsr = io.Parser('ions/', ['C_SZ_6.5au', 'C_DZDP_6.5_2.5au', 'C_TZTP_6.5_4.5_2.5au', 'H_SZ_6.5au', 'H_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au'], 'conquest/', ['CH4_SZ', 'C6H6_SZ', 'C6H6_DZDP', 'C6H6_TZTP'])
prsr.parseIons()
prsr.parseConq()
atoms = prsr.atoms
fl = prsr.fermiLevels
electrons = prsr.electrons

molecule = 'C6H6_SZ'

C = cell.Cell(molecule, fl[molecule], electrons[molecule], 25.0, 20.0, 25.0, gridSpacing=.5)
for a in atoms[molecule]:
	C.addAtom(atoms[molecule][a], a)

for i in range(0, len(C.bands)):
	E = C.bands[i]
	if abs(E - 0.0) < 0.02:
		N = i
		break
print C.electrons
C.normaliseBand(N, debug=True)
io.plot.plotChargeDensity3D(C, N, xrange=(0.0, 20.0), zrange=(0.0, 20.0), fraction=0.06, alpha=0.5, cmap=True)

#pltr.plotBasis3D('C_SZP_6.5au', 1, 3, 2, 0, show=False)
#plt.show()

# pltr.plotSPH3D(2, 0, offset=0.0)
