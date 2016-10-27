import math

# Location and names of .ion files
ionfolder = 'ions/'
ionfiles = ['H_SZ_6.5au',	'H_SZP_6.5au', 'H_DZDP_6.5_2.5au', 'H_TZTP_6.5_4.5_2.5au',
			'C_SZ_6.5au',	'C_SZP_6.5au', 'C_DZDP_6.5_2.5au', 'C_TZTP_6.5_4.5_2.5au',
			'Si_SZ_8bohr',	'Si_TZ_8_6_4bohr']

Radials = {} # Stores all data on PAO

for ion in ionfiles:
	# Open file and initialise entry
	f = open(ionfolder+ion+'.ion', 'r')
	Radials[ion] = {}

	# Skip preamble and first 9 lines
	line = f.next()
	while not '</preamble>' in line:
		line = f.next()
	for i in range(0, 9):
		line = f.next()

	# Parse PAO data
	line = f.next()
	while line.split()[0] != '#':
		metadata = line.split()
		l = metadata[0]
		n = metadata[1]
		z = metadata[2]

		orbitaldata = z+n+l
		line = f.next()
		metadata = line.split()
		pts = int(metadata[0])
		cutoff = float(metadata[2])

		Radials[ion][orbitaldata] = [[], [], cutoff]

		for i in range(0, pts):
			line = f.next()
			a, b = line.split()
			Radials[ion][orbitaldata][0].append(float(a))
			Radials[ion][orbitaldata][1].append(float(b)*math.pow(float(a),int(l)))

		line = f.next()
	f.close()