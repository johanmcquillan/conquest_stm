
import packages.plot as plot
from packages.auto_parser import Parser
import numpy as np

prsr = Parser()
C = prsr.make_cell('Si_SZP_K4_2', grid_spacing=0.5, debug=True)
del prsr

fl = C.fermi_level

# plot.plot_ldos_3d(C, -2, 0, 0.005, fraction=0.025, debug=True, top_down=False, vectorised=True)

# delta_s = 2 / (C.grid_spacing * C.H_BAR) * np.sqrt(4.5 / (2 * C.ELECTRON_MASS))
# delta_s = 15

# C.current_scan(35, -2, 0.005, 4.85, C.fermi_level, delta_s, debug=True)

z = [31, 31.5, 32, 32.5, 33, 33.5, 34, 34.5, 35]
V = np.linspace(-4.0, 4.0, 41)
for i in range(len(V)):
	for j in range(len(z)):
		# c = C.get_current_scan_plane(37, V[i], 0.005, 4.85, C.fermi_level, 32, recalculate=False, debug=True)
		plot.plot_current_2d_plane(C, 37, z[j], V[i], 0.005, 4.85, fl, 0.0, show=False, debug=True)

# G = C.greens_function_mesh(70, 4.85, C.fermi_level, debug=True)

# fraction = 0.01
#
# ld1 = C.get_ldos_grid(fl-2, fl, 0.005, debug=True)

# ld = np.ma.masked_array(ld1, mask=np.where(ld1 == 0, 1, 0))
#
# isovalue = np.max(ld)*fraction
#
#S = np.empty_like(ld)
#S[~ld.mask] = np.log(ld[~ld.mask] / isovalue)
#
#
# K = C.bands.keys()[0]
# E = C.bands[K][0]
#
# wf = C.get_psi_grid(K, E, debug=True)
# c = C.get_c(ld1, False, fraction, delta_s)
# c_mag = C.mesh_dot_product(c, c)
# print np.max(c_mag)
# print type(c_mag)
#
# #print np.min(c_mag)
#
# plot.plot_3d_scatter_mask('Test', c_mag)

