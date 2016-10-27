from parser import Parser
from plotter import Plotter
from ion import *

ionfolder = 'ions/'
ionfiles = ['H_SZ_6.5au',	'H_SZP_6.5au', 'H_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au',
			'C_SZ_6.5au',	'C_SZP_6.5au', 'C_DZDP_6.5_2.5au', 'C_TZTP_6.5_4.5_2.5au',
			'Si_SZ_8bohr',	'Si_TZ_8_6_4bohr']

# Initialise parser and get data
Prsr = Parser(ionfolder, ionfiles)
Prsr.parse()

# Put Ion objects into a dict indexed by name in ionfiles
ions = {}
for ion in ionfiles:
	ions[ion] = Prsr.getIon(ion)

# Plot radial functions to 'Rnl_radials.pdf'
Pltr = Plotter('Rnl', ions)
Pltr.plotRadials()
