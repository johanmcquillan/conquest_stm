from parser import Parser
from plotter import Plotter
from ion import *
from vector import *

n = Vector(0,0,1)
#Ion.plotSPH(1, 1, n)

ionfolder = 'ions/'
ionfiles = ['H_SZ_6.5au',	'H_SZP_6.5au', 'H_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au',
			'C_SZ_6.5au',	'C_SZP_6.5au', 'C_DZDP_6.5_2.5au', 'C_TZTP_6.5_4.5_2.5au',
			'Si_SZ_8bohr',	'Si_TZ_8_6_4bohr']

# Initialise parser and get data
Prsr = Parser(ionfolder, ionfiles)
Prsr.parseIons()

# Put Ion objects into a dict indexed by name in ionfiles
ions = {}
for ion in ionfiles:
	ions[ion] = Prsr.getIon(ion)

I = ions['C_TZTP_6.5_4.5_2.5au']
#I.plotSPH(2, 1)
I.plotBasis(1, 3, 2, 0)

# # Plot radial functions to 'Rnl_radials.pdf'
# Pltr = Plotter('Rnl', ions)
# Pltr.plotRadials()

