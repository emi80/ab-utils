#!/usr/bin/env python

import sys
import re
import itertools

inF = sys.stdin

	
def parseExons(starts, ends, coord, conn, strand):
	if not bool(conn):
		return
	if strand == "+":
		for i,c in enumerate(conn):
			if i == 0 and c == "^":
				starts.append(l[3])
				ends.append(coord[i])
				continue
			if c == "[":
				starts.append(coord[i])
			if c == "^":
				ends.append(coord[i])
			if c == "-":
				starts.append(coord[i])
			if c == "]":
				if i == 0:
					starts.append(l[3])
				ends.append(coord[i])
		if c == "-" or c == "[":
			ends.append(l[4])
		return
	if strand == "-":
		for i,c in enumerate(conn):
			if i == 0 and c == "^":
				starts.append(coord[i])
				ends.append(l[4])
				continue
			if c == "[":
				ends.append(coord[i])
			if c == "^":
				starts.append(coord[i])
			if c == "-":
				ends.append(coord[i])
			if c == "]":
				if i == 0:
					ends.append(l[4])
				starts.append(coord[i])
		if c == "-" or c == "[":
			starts.append(l[3])
		return
		

for line in inF:
	l = line.strip().strip(";").split("\t")
	d = dict(el.split(" ") for el in l[8].split("; "))

	chr, ev, strand = l[0], l[2], l[6]
	
	starts, ends = [], []

	splices = d["splice_chain"].strip("\"").split(",")
	for splice in splices:
		coord = filter(None, re.split("[\[\]\^\-]", splice))
		conn = filter(None, re.split("[0-9]+", splice))
		parseExons(starts, ends, coord, conn, strand)
	# Add encompassing exon for intron retention events
	if not splices[0]:
		conn = filter(None, re.split("[0-9]+", splices[1]))
		if len(conn) %2 == 0 and list(set(zip(*2*[iter(conn)]))) == [("^", "-")]:
			starts.append(l[3])
			ends.append(l[4])

	# Print output
	outKeys = "structure,event_id"
	outValues = "\t".join(d[key].strip("\"") for key in outKeys.split(","))
	# The coordinates are in gtf format (offset 1)
	#for s,e in itertools.izip(starts1+starts2, ends1+ends2):
	for s,e in itertools.izip(starts, ends):
		print "\t".join([chr, s, e, strand, ev, outValues])
	
		


