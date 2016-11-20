
import packages.plot as plot
from packages.parser import Parser
from packages.cell import Cell

molecule = 'CH4_TZTP'

prsr = Parser('ions/', ['Si_SZ_8bohr', 'C_SZ_6.5au', 'C_SZP_6.5au', 'C_DZDP_6.5_2.5au', 'C_TZTP_6.5_4.5_2.5au', 'H_SZ_6.5au', 'H_SZP_6.5au', 'H_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au'], 'conquest/', [molecule])
prsr.parseIons()
prsr.parseConquestOutput()
atoms = prsr.atoms
fl = prsr.fermiLevels
electrons = prsr.electrons
del prsr


C = Cell(molecule+'_161120A', fl[molecule], electrons[molecule], 15.0, 15.0, 15.0, gridSpacing=.5)
for a in atoms[molecule]:
	C.addAtom(atoms[molecule][a], a)

#C.combineBands((0.05, 0.1), normalise=True, debug=True)

for i in range(0, len(C.bands)):
	plot.plotChargeDensity3D(C, i, fraction=0.1, cmap=True, save=True, show=False)

#pltr.plotBasis3D('C_SZP_6.5au', 1, 3, 2, 0, show=False)
#plt.show()

# pltr.plotSPH3D(2, 0, offset=0.0)
