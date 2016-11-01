
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ..atomic import *

class Plotter:

	"""Stores dict of Ion objects and provides methods to plot data to pdf"""

	spectral = {0 : 's',
				1 : 'p',
				2 : 'd',
				3 : 'f'}

	def __init__(self, filename, ions):
		self.fname = filename
		self.ions = ions

	def plotRadials(self):
		"""Plot all radial functions from self.ions to 'self.filename'_radials.pdf"""
		with PdfPages('pdfs/'+self.fname+'_radials.pdf') as pdf:

			# Plot all functions for the same ion on one graph
			ionNames = sorted(self.ions.keys())
			for ionName in ionNames:
				ion = self.ions[ionName]

				# Setup plot
				plt.title('PAOs for '+ionName+'.ion')
				plt.xlabel('Radial Distance, $r$ / $a_0$')
				plt.ylabel('$R_{nl}(r)$')
				plt.minorticks_on()
				plt.grid(b=True, which='major', alpha=0.45, linestyle='-')
				plt.grid(b=True, which='minor', alpha=0.10, linestyle='-')

				# Loop over all radial functions
				for z in range(1, ion.zetas+1):
					for n in ion.nl.keys():
						for l in ion.nl[n]:
							# Get Radial and data from ion
							radial = ion.getRadial(z, n, l)
							r = radial.r
							R = radial.R
							# Add radial info to legend and add to plot
							label = '$\zeta ='+str(z)+'$, $n='+str(n)+'$, $l='+str(l)+'$'
							label = '$\zeta ='+str(z)+'$, $'+str(n)+self.spectral[l]+'$'
							plt.plot(r, R, label=label)
				# Add plot to pdf and reset plt
				plt.legend()
				pdf.savefig()
				plt.close()

	








