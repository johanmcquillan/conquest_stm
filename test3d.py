
from packages import io, cell
import matplotlib.pyplot as plt

prsr = io.Parser('ions/', ['C_SZ_6.5au', 'C_DZDP_6.5_2.5au', 'C_TZTP_6.5_4.5_2.5au', 'H_SZ_6.5au', 'H_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au'], 'conquest/', ['CH4_SZ', 'C6H6_SZ', 'C6H6_DZDP', 'C6H6_TZTP'])
prsr.parseIons()
prsr.parseConq()
atoms = prsr.atoms
fl = prsr.fermiLevels

molecule = 'CH4_SZ'

C = cell.Cell(molecule, fl[molecule], 15.0, 15.0, 15.0, gridSpacing=0.5)
for a in atoms[molecule]:
	C.addAtom(atoms[molecule][a], a)

for i in range(0, len(C.bands)):
	E = C.bands[i]
	if abs(E - 0.0) < 0.3:
		N = i
		break

io.plot.plotChargeDensity3D(C, N, fraction=0.1, alpha=0.5, cmap=True)

#pltr.plotBasis3D('C_SZP_6.5au', 1, 3, 2, 0, show=False)
#plt.show()

# pltr.plotSPH3D(2, 0, offset=0.0)
