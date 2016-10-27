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