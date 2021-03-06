#!/usr/bin/env python2.7

import conquest_stm.plot as plot
from conquest_stm.parser import get_cell

__author__ = 'Johan G. McQuillan'
__email__ = 'johan.mcquillan.13@ucl.ac.uk'

debug = True

C = get_cell('Si001_test', grid_spacing=1.0, debug=debug)

fl = C.fermi_level
V = -2
T = 10
z = 40.0

plot.plot_ldos_3d(C, V, 0, T, fraction=0.01, debug=debug)
plot.plot_current_2d(C, z, V, T, 4.85, fl, 0.01, debug=debug)
