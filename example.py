
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
for ion in ionFilesRaw:
	ions[ion] = Prsr.getIon(ion)

# ions is a dict of all the Ion objects, indexed by the names given in ionFilesRaw
# An Ion object 











