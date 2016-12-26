
import packages.plot as plot
from packages.auto_parser import Parser

prsr = Parser()
C = prsr.make_cell('Si_SZP_K4', grid_spacing=0.5, debug=True)
del prsr

#C.write_ldos(-2, 0, 0.01, debug=True)

#plot.plot_charge_density_gamma_3d(C, 0)

plot.plot_ldos_3d(C, -2, 0, 0.005, fraction=0.01, save=False, show=True, debug=True)
#	for z in [29.5, 30, 30.5, 31, 31.5, 32]:
#		plot.plot_ldos_2d(C, -2.0, 0.0, 0.005, 'z', z, 0, 14, bias=bias, debug=True)