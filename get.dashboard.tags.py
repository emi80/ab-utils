#!/usr/bin/env python


import sys



if len(sys.argv) != 3:
	print "USAGE: %s index_file.txt <tags>" %sys.argv[0]
	print "<tags> should be valid gtf keys in the 9th field, comma-separated"
	exit()


input = sys.stdin if sys.argv[1] == "-" else open(sys.argv[1])


which = sys.argv[2].split(",")
for line in input:
	f,m = line.strip().strip(";").split("\t")
	d = dict(el.split("=") for el in m.split("; "))
	l = []
	for w in which:
		if w=="path":
			l.append(f)
			continue
		l.append(d.get(w, "NA"))
	print "\t".join(l)
input.close()
exit()
