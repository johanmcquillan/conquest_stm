from parser import Parser
from plotter import Plotter
from ion import *
from vector import *

ionfolder = 'ions/'
# ionfiles = ['H_SZ_6.5au',	'H_SZP_6.5au', 'H_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au',
# 			'C_SZ_6.5au',	'C_SZP_6.5au', 'C_DZDP_6.5_2.5au', 'C_TZTP_6.5_4.5_2.5au',
# 			'Si_SZ_8bohr',	'Si_TZ_8_6_4bohr']



# ionfiles = ['C_TZTP_6.5_4.5_2.5au']
# Initialise parser and get data
Prsr = Parser(ionfolder, ionfiles)
Prsr.parseIons()

# Put Ion objects into a dict indexed by name in ionfiles
ions = {}
for ion in ionfiles:
	ions[ion] = Prsr.getIon(ion)

#Pltr = Plotter('Rnl', ions)
#Pltr.plotRadials()

I = ions['C_TZTP_6.5_4.5_2.5au']
#I.plotSPH(2, 1)
for n in I.nl.keys():
	for l in I.nl[n]:
		for m in range(-l, l+1):
			for e in ['x','y','z']:
				#I.plotBasis(1, n, l, m, e)
				continue

for l in range(0, 4):
	for m in range(-l, l+1):
		for e in ['x', 'y', 'z']:
			Ion.plotSPH(l, m, e)

# # Plot radial functions to 'Rnl_radials.pdf'
# Pltr = Plotter('Rnl', ions)
# Pltr.plotRadials()

