
from packages.parser import Parser
import packages.plot as plot

ionFolder = 'ions/'
ionFiles = ['H_SZ_6.5au', 'C_SZ_6.5au', 'H_SZP_6.5au', 'C_SZP_6.5au', 'H_DZDP_6.5_2.5au', 'C_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au', 'C_TZTP_6.5_4.5_2.5au']

prsr = Parser(ionFolder, ionFiles, 'conquest/', [])
prsr.parseIons()

ions = prsr.ions

plot.plotRadials(ions, spectro=True)
