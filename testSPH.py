
from packages import io

lMax = 3

print "3-DIMENSIONAL SPHERICAL HARMONICS"
print "   Enter 'q' to quit"

quit = False
while not quit:
	gotL = False
	while not gotL and not quit:
		lString = raw_input("Enter Degree, l = ")
		if lString == "q":
			quit = True
		else:
			try:
				l = int(lString)
				if l > lMax or l < 0:
					raise ValueError
				else:
					gotL = True
			except ValueError:
				print "l must be between an integer between 0 and 3 inclusive"

	if not quit:
		if l == 0:
			m = 0
		else:
			gotM = False
			while not gotM and not quit:
				if l == 0:
					m = 0
				else:
					mString = raw_input("Enter Order, m = ")
					if mString == "q":
						quit = True
					else:
						try:
							m = int(mString)
							if m > l or m < -l:
								raise ValueError
							else:
								gotM = True
						except ValueError:
							print "m must be and integer "+str(-l)+" and "+str(+l)+" inclusive"
		if not quit:
			print "Showing Spherical Harmonic for l="+str(l)+" and m="+str(m)
			io.plot.plotSPH3D(l, m)
