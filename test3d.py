
import packages.plot as plot
from packages.auto_parser import Parser
from packages.vector import KVector

molecule = 'Si_SZP_K4'
molecule = 'C6H6_DZDP'
prsr = Parser()

ion = prsr.get_ion('C_DZDP_6.5_2.5au')
plot.plot_radials({ion.ion_name: ion})
# C = prsr.make_cell('C6H6_DZDP', grid_spacing=0.25, debug=True)
# fl = C.fermiLevel
# plot.plot_ldos_3d(C, -2, 0, 0.05, debug=True)

# for i in range(0, len(C.bands)):
	# plot.plotChargeDensity3D(C, i, fraction=0.1, cmap=True, save=True, show=False)

#plot.plotChargeDensity3D(cell, cell.fermiLevel, cell.fermiLevel-1.0, fraction=0.05, cmap=True, save=False, show=True, debug=True)

#plot.plot_ldos_2d(cell, -2, 0, 0.01, 'z', 0, 14, 32, interpolation='quadratic', debug=True)

