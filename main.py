
import sys

import packages.plot as plot

from argparse import ArgumentParser
from packages.auto_parser import Parser

arg_parser = ArgumentParser(description='Get input filename')
arg_parser.add_argument('input', type=str, help='Input file to read')

input_filename = vars(arg_parser.parse_args())['input']
input_file = open(input_filename, 'r')


settings = {}

end_of_file = False
while not end_of_file:
	try:
		line = input_file.next()
		line_partition = line.partition('#')  # Partition by comment symbol
		line_split = line_partition[0].split()  # Get everything before comment symbol
		if len(line_split) > 1:
			setting = line_split[0]
			argument = tuple(line_split[1:])
			settings[setting] = argument
	except StopIteration:
		end_of_file = True


plot_methods = {'plot_sph_3d', 'plot_ldos_2d', 'plot_ldos_3d'}

if 'Conquest_out' not in settings:
	print 'Need Conquest_out field'
	sys.exit(1)
else:
	prsr = Parser()
	cell = prsr.make_cell(str(settings['Conquest_out'][0]))

	for s in settings:
		if s != 'Conquest_out':
			print s[5:]
			if s.startswith('plot.') and s[5:] in plot.__dict__:
				plot.__dict__[s[5:]](*settings[s])
