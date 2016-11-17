
# from packages import cell
# from packages import io
# from packages import atomic

from packages.parser import Parser
from packages.cell import Cell
import packages.plot as plot

# Folder of the .ion files.
ionFolder = 'ions/'
# ion files to be read, excluding .ion
ionFiles = ['H_SZ_6.5au', 'C_SZ_6.5au']

# Folder of Conquest output files
conqFolder = 'conquest/'
# Name of simulation files
conqFile = 'C6H6_SZ'
# This denotes the files C6H6_SZ (Conquest_out), C6H6_SZ.dat (Process0000001WF.dat),
#  and C6H6_SZ.dos (DOSout.dat)

# To parse the data, instantiate a Parser object.
prsr = Parser(ionFolder, ionFiles, conqFolder, [conqFile])

# First, parse ions into Ion objects.
#  These hold all of the radial functions, indexed by zeta, n, and l
prsr.parseIons()
# Then combine with Conquest output
prsr.parseConquestOutput()
# Parser now holds Atom objects, which extend the Ion class and stores
#  atomic position and coefficients

# Currently, atoms, Fermi levels, and total number of electrons must be taken
#  from Parser and put into a cell, but I plan to write a method to get the cell
#  directly from Parser

# Parser can hold data from several simulations at once, so when getting data
#  from it, they are indexed by simulation name
atoms = prsr.atoms[conqFile] # Dict of Atom objects
fermi = prsr.fermiLevels[conqFile]
electrons = prsr.electrons[conqFile]

# Create a 3D cell with dimensions 20 x 20 x 15 a0
cell = Cell(conqFile+'_EXAMPLE', fermi, electrons, 20.0, 20.0, 15.0, gridSpacing=.5)

# Add the atoms to the cell
for atomKey in atoms:
	cell.addAtom(atoms[atomKey], atomKey)

# io.plot has methods to plot 2D and 3D

# For 7th band, plot charge density isosurface with surface value of 5% of max value
plot.plotChargeDensity3D(cell, 6, fraction=0.05, show=True, save=True, cmap=True)

# C.bands is an ordered list of band energies for the cell
#  This loop, when activated, will save 3D plots of all bands to folder figures3D
#  (You may have to create the folder locally first)
if False:
	for i in range(len(cell.bands)):
		plot.plotChargeDensity3D(cell, i, fraction=0.05, show=False, save=True, cmap=True)
