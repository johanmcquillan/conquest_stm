
import packages.plot as plot
from packages.auto_parser import Parser

prsr = Parser()
C = prsr.make_cell('Si_SZP_K4', grid_spacing=0.5, debug=True)
del prsr



#C.write_ldos_grid(-2, 0, 0.01, debug=True)

#plot.plot_charge_density_gamma_3d(C, 0)

fl = C.fermi_level
print fl
print C.bands[C.bands.keys()[0]][-1]
plot.plot_ldos_3d(C, -2, 0, 0.005, fraction=0.05, debug=True, top_down=False, vectorised=True)

#C.calculate_support_grid(vectorised=True, debug=True)

#K = C.bands.keys()[0]
#E = C.bands[K][50]
#C.calculate_psi_grid(K, E, recalculate=False, write=False, debug=True, vectorised=False)
# C.calculate_psi_grid(K, E, recalculate=False, write=False, debug=True, vectorised=True)

#C.bardeen_element(C.c(C.get_psi_range(-2, 0, 0.005, debug=True), 0.0005, 0.001), )
#plot.plot_vector_field(C, C.c(C.get_ldos_grid(-2, 0, 0.005, debug=True), 0.0005, 0.001))

#C.broadened_volume_integrand(C.get_ldos_grid(-2+fl, fl, 0.005, debug=True), 1, 0.05)

#plot.plot_ldos_2d(C, -2, 0, 0.005, 'z', 30.5, 0, 14, debug=True)

#plot.plot_ldos_3d(C, -2, 0, 0.005, fraction=0.05, save=False, show=True, debug=True)
#	for z in [29.5, 30, 30.5, 31, 31.5, 32]:
#		plot.plot_ldos_2d(C, -2.0, 0.0, 0.005, 'z', z, 0, 14, bias=bias, debug=True)
