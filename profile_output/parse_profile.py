
profile = open("profile_161130_160647.dat")

profile.next()
profile.next()
profile.next()
profile.next()
profile.next()

calls = []

for line in profile:
	calls.append(line.split())

calls = sorted(calls, key=lambda x: float(x[3]), reverse=True)

for i in range(20):
	toPrint = ""
	for j in range(len(calls[i])):
		toPrint += calls[i][j]
		if j != len(calls[i]) - 1:
			toPrint += "   "
	print toPrint
