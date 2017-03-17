#!/usr/bin/env python

import argparse
import sys
import itertools

parser = argparse.ArgumentParser(description='Annotate the membership of elements to a given set')
parser.add_argument('-i', '--input', type=str, default="stdin", help='Two-column file: element, set. [default=stdin]')
#parser.add_argument('-u', '--up', type=int, default='0', help='How many nucleotides to extend upstream')
parser.add_argument('-p', '--pseudo', type=str, 
	help='One or more labels comma-separated for sets. If the label is already in the input it is ignored, otherwise it will be an empty set')
parser.add_argument('--combs', action='store_true', default=False, help='Output the membership to all possible combinations of sets')
parser.add_argument('--bool', action='store_true', default=False, help="Write true or false instead of 0 and 1")
args = parser.parse_args()


input = sys.stdin if args.input=="stdin" else open(args.input)

bool_dict = {"0": "0", "1": "1"}
if args.bool:
	bool_dict = {"0": "false", "1": "true"}

d = {}
sets = set()
for line in input:
	l = line.strip().split("\t")
	d.setdefault(l[0], {})
	d[l[0]][l[1]] = 1
	sets.add(l[1])

if args.pseudo:
	sets.update(args.pseudo.split(","))

header = "element\t" + "\t".join(sorted(sets)) + "\tsum" + "\tmembership"

combs = ""
if args.combs:
	combs = sum(list( list(itertools.combinations(sets, r)) for r in range(2, len(sets)+1) ), [])
	header = header + "\t" + "\t".join( list(",".join(comb) for comb in combs ) )
print header

for k,v in sorted(d.iteritems()):
	binary = "\t".join( bool_dict[str(v.get(s, 0))] for s in sorted(sets))
	nb_sets = str(sum(v.itervalues()))
	memb = ",".join(sorted(v.keys()))
	line = "\t".join( [ k, binary, nb_sets, memb ] )
	if args.combs:
		line_combs = "\t".join( bool_dict[str ( int( set(comb) == set(v) ) )] for comb in combs )
		line = line + "\t" + line_combs
	print line
			
exit()
