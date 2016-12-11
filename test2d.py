
import packages.plot as plot
from packages.auto_parser import Parser

prsr = Parser()
C = prsr.make_cell('Si_SZP_K4', grid_spacing=0.5, debug=True)
del prsr

#C.write_ldos(-2, 0, 0.01, debug=True)

plot.plot_ldos_3d(C, -2, 0, 0.005, fraction=0.05, debug=True)