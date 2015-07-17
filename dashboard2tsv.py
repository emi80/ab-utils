#!/usr/bin/python

import sys


if len(sys.argv) == 1:
	print "USAGE: %s <dashboard_index.txt> <output_fname.tsv> \n Leave '-' if you want to read the index file from stdin" %(__file__)
	exit()


input_fname = sys.argv[1]
output_fname = sys.argv[2]

open_input = sys.stdin if input_fname == "-" else open(input_fname)

d = {}

for i,line in enumerate(open_input):
	if line.startswith('#'):
		continue
	fname, tags = line.strip('\n; ').split('\t')
	tags = dict((tag.strip().split('=')[0], "=".join(tag.strip().split('=')[1:])) for tag in tags.strip().split(';'))
	d[fname] = tags

header = list(sorted(set(key for tags in d.itervalues() for key in tags.iterkeys())))
#header = list(sorted(set(key for tags in d.itervalues() for key in tags.iterkeys())))

output_file = open(output_fname, 'w')
output_file.write('\t'.join(["path"] + header)+'\n')
for key, values in d.iteritems():
	line = '\t'.join([key] + map(lambda x: values.get(x,'NA').strip("\""), header))
	output_file.write(line+'\n')
output_file.close() 

