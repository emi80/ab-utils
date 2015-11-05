#!/usr/bin/env python

# Import modules
from argparse import ArgumentParser
from subprocess import Popen, PIPE
import sys

# Argument parsing

parser = ArgumentParser(description='Parse an obo file and produces two tsv files: 1) parent-child relations 2) term description')
parser.add_argument("-i", "--input", type=str, default="stdin", help="obo file. \"stdin\" to read from stdin")
parser.add_argument("-o", "--output", type=str, default="obo.out", help="output file name WITHOUT extension")
#parser.add_argument("-g", "--gtf", type=str, help="gtf file with all elements (genes, transcripts, exons).")
#parser.add_argument("-w", "--window", type=int, help="number of nucleotides around junctions.", default=20)
args = parser.parse_args()

#

def read_obo(fo):
	start = False
	d = {}
	for line in fo:
		if line.startswith("#"):
			continue
		if line.startswith("[Term]"):
			start = True
			continue
		if line.strip() == "":
			start = False
			continue	
		if start:
			if line.startswith("id"):
				id = line.strip().split()[1]
				d[id] = {}
			else:
				try:
					line_split = line.strip().split(" ! ")[0]
					key = line_split.split(": ")[0]
					value = " ".join(line_split.split(": ")[1:])
				except:
					print line
					print line_split
					exit()
				d[id].setdefault(key, []).append(value)
	return d




if args.input == "stdin":
	d = read_obo(sys.stdin)
else:
	d = read_obo(open(args.input))


out1 = open(args.output + ".relationship.tsv", "w")
out2 = open(args.output + ".description.tsv", "w")


for k, v in d.iteritems():
	# Ignore obsolete terms
	if v.get("is_obsolete", 0) == ["true"]:
		continue
	out2.write("%s\t%s\t%s\n" %(k, v["name"][0], v["namespace"][0]))
	for rel in v.get("is_a", ["NA"]):
		out1.write("%s\t%s\n" %(rel, k))

out1.close()
out2.close()

