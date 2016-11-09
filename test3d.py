
from packages import io
from packages import cell
import matplotlib.pyplot as plt

prsr = io.Parser('ions/', ['C_SZ_6.5au', 'H_SZ_6.5au'], 'conquest/', ['CH4_SZ', 'C6H6_SZ'])
prsr.parseIons()
prsr.parseConq()
atoms = prsr.atoms

pltr = io.Plotter('test', prsr.ions)


C = cell.Cell('CH4_SZ', 15.0, 15.0, 15.0, gridSpacing=0.2)
for a in atoms['CH4_SZ']:
	C.addAtom(atoms['CH4_SZ'][a], a)


pltr.plotChargeDensity3D(C, 0, minimum=-15.0, maximum=+15.0)

#pltr.plotBasis3D('C_SZP_6.5au', 1, 3, 2, 0, show=False)
#plt.show()

# pltr.plotSPH3D(2, 0, offset=0.0)
