import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class Plotter:

	def __init__(self, filename, ions):
		self.fname = filename
		self.ions = ions

	def plotRadials(self):
		# Plot functions to pdf
		with PdfPages(self.fname+'_radials.pdf') as pdf:
			ionNames = sorted(self.ions.keys())
			for ionName in ionNames:
				ion = self.ions[ionName]
				for z in range(1, ion.zetas+1):
					for n in ion.nl.keys():
						print ionName, ion.nl[n]
						for l in ion.nl[n]:
							label = '$\zeta ='+str(z)+'$, $n='+str(n)+'$, $l='+str(l)+'$'
							radial = ion.getRadial(z, n, l)
							r = radial.r
							R = radial.R
							plt.plot(r, R, label=label)
				plt.title('PAOs for '+ionName+'.ion')
				plt.xlabel('Radial Distance, $r$ / $a_0$')
				plt.ylabel('$R_{nl}(r)$')
				plt.minorticks_on()
				plt.grid(b=True, which='major', alpha=0.45, linestyle='-')
				plt.grid(b=True, which='minor', alpha=0.10, linestyle='-')
				plt.legend()
				pdf.savefig()
				plt.close()
