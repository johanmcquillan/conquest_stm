import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math

# Folder and names
ionfolder = 'ions/'
ionfiles = ['H_SZ_6.5au',	'H_SZP_6.5au', 'H_DZDP_6.5_2.5au',
			'C_SZ_6.5au',	'C_SZP_6.5au', 'C_DZDP_6.5_2.5au',
			'Si_SZ_8bohr']

Radials = {}

for ion in ionfiles:
	f = open(ionfolder+ion+'.ion', 'r')

	Radials[ion] = {}

	line = f.readline()

	while not '</preamble>' in line:
		line = f.readline()

	for i in range(0, 9):
		line = f.readline()

	line = f.readline()

	while line.split()[0] != '#':
		print ion, line
		metadata = line.split()
		l = metadata[0]
		n = metadata[1]
		z = metadata[2]

		orbitaldata = z+n+l
		print orbitaldata
		line = f.readline()
		metadata = line.split()
		pts = int(metadata[0])
		cutoff = float(metadata[2])

		Radials[ion][orbitaldata] = [[], [], cutoff]

		for i in range(0, pts):
			line = f.readline()
			a, b = line.split()
			Radials[ion][orbitaldata][0].append(float(a))
			Radials[ion][orbitaldata][1].append(float(b))

		line = f.readline()

with PdfPages('Rnl.pdf') as pdf:
	for ion in ionfiles:
		orbitals = Radials[ion].keys()
		orbitals = sorted(orbitals, key=lambda x: x)
		cutoffs = []
		for o in orbitals:
			cutoffs.append(Radials[ion][o][2])
			label = '$\zeta ='+o[0]+'$, $n='+o[1]+'$, $l='+o[2]+'$'#, $m_l='+o[2]+'$'
			plt.plot(Radials[ion][o][0], Radials[ion][o][1], label=label)
		plt.xlim(0, math.ceil(max(cutoffs)))
		plt.title('Radial Functions for '+ion+'.ion')
		plt.xlabel('Radius, $r$ / $a_0$')
		plt.ylabel('$R_{nl}(r)$')
		plt.grid(b=True, which='both', color='0.65',linestyle='--')
		plt.legend()
		pdf.savefig()
		plt.close()