
import packages.plot as plot
from packages.parser import Parser

prsr = Parser('ions/', ['Si_SZ_8bohr'], 'conquest/', ['Si_SZ'])
prsr.parseIons()
prsr.parseConquestOutput()
C = prsr.getCell('Si_SZ')

C.name+= '_161119A'

plot.plotChargeDensity2D(C, 5, 'z', minimum=-10.0, maximum=10.0)
