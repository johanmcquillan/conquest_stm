import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
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

# Plot functions to pdf
with PdfPages('Rnl.pdf') as pdf:
	for ion in ionfiles:
		orbitals = sorted(Radials[ion].keys())
		cutoffs = []
		for o in orbitals:
			cutoffs.append(Radials[ion][o][2])
			label = '$\zeta ='+o[0]+'$, $n='+o[1]+'$, $l='+o[2]+'$'#, $m_l='+o[2]+'$'
			plt.plot(Radials[ion][o][0], Radials[ion][o][1], label=label)
		plt.xlim(0, math.ceil(max(cutoffs)))
		plt.title('PAOs for '+ion+'.ion')
		plt.xlabel('Radial Distance, $r$ / $a_0$')
		plt.ylabel('$R_{nl}(r)$')
		plt.minorticks_on()
		plt.grid(b=True, which='major', alpha=0.45, linestyle='-')
		plt.grid(b=True, which='minor', alpha=0.10, linestyle='-')
		plt.legend()
		pdf.savefig()
		plt.close()