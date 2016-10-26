import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


ionfiles = ['H_SZ_6.5au', 'C_SZ_6.5au']

x = {}
y = {}

for ion in ionfiles:
	f = open(ion+'.ion', 'r')

	foundPAO = False

	line = f.readline()

	while not '</preamble>' in line:
		line = f.readline()

	for i in range(0, 9):
		line = f.readline()

	metadata = f.readline().split()
	l = metadata[1]
	n = metadata[2]
	z = metadata[3]

	metadata = f.readline().split()
	pts = int(metadata[0])

	x[ion] = []
	y[ion] = []

	for i in range(0, pts):
		line = f.readline()
		a, b = line.split()
		x[ion].append(float(a))
		y[ion].append(float(b))

with PdfPages('Rnl.pdf') as pdf:
	for ion in ionfiles:
			plt.plot(x[ion], y[ion])
			plt.title('Radial Part for '+ion+'.ion')
			plt.xlabel('Radius, $r$ / $a_0$')
			plt.ylabel('$R_{nl}(r)$')
			plt.grid(b=True, which='both', color='0.65',linestyle='--')
			pdf.savefig()
			plt.close()