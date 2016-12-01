
import packages.plot as plot
from packages.parser import Parser

prsr = Parser('ions/', ['C_SZ_6.5au', 'H_SZ_6.5au'], 'conquest/', ['C6H6_SZ'])
prsr.parseIons()
prsr.parseConquestOutput()
C = prsr.getCell('C6H6_SZ', gridSpacing=0.5)

plot.plotLDoS2D(C, -2, 0, 0.01, 'z', 0, 20, interpolation='cubic', debug=True)
