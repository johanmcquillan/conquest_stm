
import glob
import packages.io as IO
import packages.atomic as atm

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
Prsr = IO.Parser()

# Currently, Parser only works for .ion files
Prsr.parseIons(ionFolder, ionFilesAll)

# Prsr now holds all the info about the ions.
# In the final code, we would run something like
#  Prsr.parseAtoms(atomFolder, atomFiles)
#  which would combine the ion data with the other input files to get the atoms

# For now, we get the ions from Prsr
ions = {}
for ion in ionFilesAll:
	ions[ion] = Prsr.getIon(ion)

# ions is a dict of all the Ion objects, indexed by the names given in ionFilesAll
# An Ion object stores the corresponding Radial objects, which are the PAO's given
#  in the .ion file.
# To see the radial functions, we use Plotter
Pltr = IO.Plotter('example', ions)
Pltr.plotRadials()

# Running this code will produce example_radials.pdf in the folder pdfs

# The Ion class combines the Radial data with the spherical harmonics
#  to form the basis functions. We can plot a 2D cross section of the 
#  basis functions. Currently, this is done as a method in Ion class,
#  but will be moved to the Plotter class such that output can be put
#  into a single pdf.

I = ions['C_TZTP_6.5_4.5_2.5au']

# I.nl is a dict that stores all of the l values for a particular n
# eg. I.nl[2] = [0, 1]; I.nl[3] = [2]
for n in I.nl.keys():
	for l in I.nl[n]:
		for m in range(-l, l+1):
			# e will be the axis that is set to a constant, planeValue, to get a 2D plot
			for e in ['x','y','z']:
				I.plotBasis(1, n, l, m, e, step=0.1, planeValue=0.0)

# The plots will be output to pdfs folder

# The next step will involve an Atom class that inherits from the Ion
#  class and also stores the coefficients such that the wavefunction
#  around the atom can be found.













