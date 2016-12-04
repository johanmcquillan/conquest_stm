
import argparse

parser = argparse.ArgumentParser(description='Get filename')
parser.add_argument('filename', metavar='i', type=str, help='File to read')

filename = vars(parser.parse_args())['filename']
profile = open(filename)

profile.next()
profile.next()
profile.next()
profile.next()
headers = profile.next().split()

calls = []

for line in profile:
	if len(line.split()) == 6:
		calls.append(line.split())
	else:
		calls[-1].extend(line.split())

calls = sorted(calls, key=lambda x: float(x[3]), reverse=True)

width = '13'
rowformat = '%'+width+'s %'+width+'s %'+width+'s %'+width+'s %'+width+'s %'+width+'s'

print rowformat % (headers[0], headers[1], headers[2], headers[3], headers[4], headers[5])
for i in range(20):
	toPrint = rowformat % (calls[i][0], calls[i][1], calls[i][2], calls[i][3], calls[i][4], calls[i][5])
	print toPrint
