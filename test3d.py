
import packages.plot as plot
from packages.parser import Parser
from packages.cell import Cell

molecule = 'C6H6_SZ'

prsr = Parser('ions/', ['Si_SZ_8bohr', 'C_SZ_6.5au', 'C_SZP_6.5au', 'C_DZDP_6.5_2.5au', 'C_TZTP_6.5_4.5_2.5au', 'H_SZ_6.5au', 'H_SZP_6.5au', 'H_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au'], 'conquest/', [molecule])
prsr.parseIons()
prsr.parseConquestOutput()

cell = prsr.getCell(molecule, gridSpacing=0.5)
del prsr
cell.name += '_161125'

# for i in range(0, len(C.bands)):
	# plot.plotChargeDensity3D(C, i, fraction=0.1, cmap=True, save=True, show=False)

plot.plotChargeDensity3D(cell, cell.fermiLevel, cell.fermiLevel-0.5, fraction=0.05, cmap=True, save=False, show=True)
