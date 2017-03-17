#! /usr/bin/env python

import argparse, sys

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Sort a matrix based on a vector')
parser.add_argument("-i", "--matrix", type=str, default="stdin", help="Matrix you want to sort. Can be stdin")
parser.add_argument("-c", "--colOrder", type=str, help="Ordered list of column names")
#parser.add_argument("-x", type=str, default="1", help="index of the field you want file a to be joined on. Accept multiple indeces comma-separated")
#parser.add_argument("-y", type=str, default="1", help="index of the field you want file b to be joined on. Accept multiple indeces comma-separated")
#parser.add_argument("--a_header", action="store_true", help="file \"a\" has a header")
#parser.add_argument("--b_header", action="store_true", help="file \"b\" has a header")
args = parser.parse_args()


# ===========
# BEGIN
# ===========

colOrder = [line.strip() for line in open(args.colOrder)]

matrix = sys.stdin if args.matrix == "stdin" else open(args.matrix)

sortCols = True

SEP = "\t"

for i,line in enumerate(matrix):
	if i == 0:
		header = line.strip().split(SEP)
		if sortCols:
			# Retain only ids that are in the header
			colOrder = list(el for el in colOrder if el in header)
			d = dict((colOrder.index(h.strip()),i) for i,h in enumerate(header))
			print SEP.join(header[d.get(i)] for i in range(len(header)))
		continue
	shifted = len(line.split(SEP)) == len(header)+1 
	offset = int(shifted)
	values = line.strip().split(SEP)
	if sortCols:
		print values[0] + SEP + SEP.join(values[d.get(i)+offset] for i in range(len(header)))
	
