
import packages.plot as plot
from packages.parser import Parser

prsr = Parser('ions/', ['Si_SZP_8bohr'], 'conquest/', ['Si_SZP_K4'])
prsr.parseIons()
prsr.parseConquestOutput()
C = prsr.getCell('Si_SZP_K4', gridSpacing=0.5)
del prsr

plot.plotLDoS2D(C, -2, 0, 0.01, 'z', 0, 14, 32, debug=True)
