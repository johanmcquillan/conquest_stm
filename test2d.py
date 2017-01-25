
import packages.plot as plot
from packages.auto_parser import Parser
import numpy as np

prsr = Parser()
C = prsr.make_cell('Si_SZP_K4', grid_spacing=0.5, debug=True)
del prsr

fl = C.fermi_level

plot.plot_ldos_3d(C, -2, 0, 0.005, fraction=0.05, debug=True, top_down=False, vectorised=True)

# delta_s = 2 / (C.grid_spacing * C.H_BAR) * np.sqrt(4.5 / (2 * C.ELECTRON_MASS))

# C.current_scan(35, -2, 0.005, 4.85, C.fermi_level, delta_s, debug=True)

# plot.plot_current_2d(C, 35, -2, 0.005, 4.85, C.fermi_level, delta_s, recalculate=False, debug=True)
