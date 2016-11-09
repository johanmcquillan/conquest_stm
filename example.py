
import glob

from packages import io
from packages import atomic

# First, we need to provide the folder of the .ion files.
ionFolder = 'ions/'

# We can either tell it which specific .ion files to read...
ionFilesSpecific = ['H_SZ_6.5au', 'C_SZ_6.5au'] # File names, excluding .ion
#  or we can use glob to get all .ion files within ionFolder.
ionFilesRaw = glob.glob(ionFolder+'*.ion')

# ionFilesRaw will store the full path eg. 'ion_parse/ions/H_SZ_6.5au.ion'
# We need to get rid of the folder info and the .ion extension.
ionFilesAll = []
tempName = ''
for ionNameRaw in ionFilesRaw:
	tempName = ionNameRaw[len(ionFolder):len(ionNameRaw)-4]
	ionFilesAll.append(tempName)

# ionFilesAll now holds the names of the ions like in ionFilesSpecific,
#  but for all ions within ionFolder.

# To parse the data, we instantiate a Parser object.
# This handles reading data from several different input files and locations,
#  which will be needed as we need to combine info from .ion, Conquest_out
#  and coefficient .dat files.
Prsr = io.Parser()

# Currently, Parser only works for .ion files
Prsr.parseIons(ionFolder, ionFilesAll)

# Prsr now holds all the info about the ions.
# In the final code, we would run something like
#  Prsr.parseAtoms(atomFolder, atomFiles)
#  which would combine the ion data with the other input files to get the atoms

# For now, we get the ions from Prsr
ions = Prsr.ions

# ions is a dict of all the Ion objects, indexed by the names given in ionFilesAll
# An Ion object stores the corresponding Radial objects, which are the PAO's given
#  in the .ion file.
# To see the radial functions, we use Plotter
Pltr = io.Plotter('example', ions)
Pltr.plotRadials(printStatus=True)

# Running this code will produce example_radials.pdf in the folder pdfs

# The Ion class combines the Radial data with the spherical harmonics
#  to form the basis functions. We can plot a 2D cross section of the
#  basis functions. Currently, this is done as a method in Ion class,
#  but will be moved to the Plotter class such that output can be put
#  into a single pdf.

ionName = 'C_TZTP_6.5_4.5_2.5au'
ion = ions[ionName]

for n in ion.radials[1]:
	for l in ion.radials[1][n]:
		for m in range(-l, l+1):
			# e will be the axis that is set to a constant, planeValue, to get a 2D plot
			for e in ['x', 'y', 'z']:
				Pltr.plotBasis(ionName, 1, n, l, m, e, step=0.1, planeValue=0.0, printStatus=True)

# The plots will be output to pdfs folder

# The next step will involve an Atom class that inherits from the Ion
#  class and also stores the coefficients such that the wavefunction
#  around the atom can be found.
