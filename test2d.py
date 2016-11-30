
import packages.plot as plot
from packages.parser import Parser

prsr = Parser('ions/', ['C_SZ_6.5au', 'H_SZ_6.5au'], 'conquest/', ['C6H6_SZ'])
prsr.parseIons()
prsr.parseConquestOutput()
C = prsr.getCell('C6H6_SZ', gridSpacing=1)

plot.plotLDoS2D(C, 5, 'z', minimum=0.0, maximum=20.0, debug=True)
