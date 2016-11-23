
from packages.parser import Parser
from packages.cell import Cell

prsr = Parser('ions/', ['H_SZ_6.5au', 'C_SZ_6.5au'], 'conquest/', ['C6H6_SZ'])
prsr.parseIons()
prsr.parseConquestOutput()

cell = prsr.getCell('C6H6_SZ', gridSpacing=0.05)

print cell.getTotalCharge(0)

