
from packages import io
from packages import cell
import matplotlib.pyplot as plt

prsr = io.Parser('ions/', ['C_SZ_6.5au', 'C_DZDP_6.5_2.5au', 'H_SZ_6.5au', 'H_DZDP_6.5_2.5au'], 'conquest/', ['CH4_SZ', 'C6H6_SZ', 'C6H6_DZDP'])
prsr.parseIons()
prsr.parseConq()
atoms = prsr.atoms

pltr = io.Plotter('test')


C = cell.Cell('CH4_SZ', 15.0, 15.0, 15.0, gridSpacing=0.5)
for a in atoms['CH4_SZ']:
	C.addAtom(atoms['CH4_SZ'][a], a)

for i in range(0, len(C.bands)):
	E = C.bands[i]
	if abs(E - 0.1) < 0.05:
		N = i

pltr.plotChargeDensity3D(C, int(len(C.bands)/2), xrange=(0.0, 15.0), yrange=(0.0,15.0), zrange=(0.0,15.0), step=C.gridSpacing, fraction=0.2, alpha=0.5)

#pltr.plotBasis3D('C_SZP_6.5au', 1, 3, 2, 0, show=False)
#plt.show()

# pltr.plotSPH3D(2, 0, offset=0.0)
