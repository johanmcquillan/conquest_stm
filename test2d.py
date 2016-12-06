
import packages.plot as plot
from packages.auto_parser import Parser

prsr = Parser()
C = prsr.make_cell('Si_SZP_K4', grid_spacing=0.5)
del prsr

plot.plot_ldos_2d(C, -2, 0, 0.01, 'z', 0, 14, 29, debug=True)
