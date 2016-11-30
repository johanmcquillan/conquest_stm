
import packages.plot as plot
from packages.parser import Parser
from packages.cell import Cell

molecule = 'Si_SZP_K4'

prsr = Parser('ions/', ['Si_SZP_8bohr'], 'conquest/', [molecule])
prsr.parseIons()
prsr.parseConquestOutput()

cell = prsr.getCell(molecule, gridSpacing=0.5)
del prsr
cell.name += '_161125'

# for i in range(0, len(C.bands)):
	# plot.plotChargeDensity3D(C, i, fraction=0.1, cmap=True, save=True, show=False)

#plot.plotChargeDensity3D(cell, cell.fermiLevel, cell.fermiLevel-1.0, fraction=0.05, cmap=True, save=False, show=True, debug=True)

plot.plotLDoS2D(cell, -2, 0, 0.01, 'z', 0, 14, 32, debug=True)

