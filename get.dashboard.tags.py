#!/usr/bin/env python


import sys

usage="\n\
	USAGE: %s index_file.txt <tags> \n\
	\n\
	<tags> should be valid keys in the 2nd field, comma-separated.\n\
	Use \"path\" if you want also the path as key\n\
	Add -header if you want to print the header as well\n\
" %sys.argv[0]

if len(sys.argv) < 3:
	print usage
	exit()

bn = sys.argv[0].split("/")[-1].replace(".py", "")

input = sys.stdin if sys.argv[1] == "-" else open(sys.argv[1])

try:
	if sys.argv[3] == "-header":
		print "\t".join(sys.argv[2].split(","))
	else:
		print usage
		exit()	
except:
	pass

which = sys.argv[2].split(",")
for line in input:
	f,m = line.strip().strip(";").split("\t")
	if m.startswith("#"):
		continue
	try: 
		d = dict(el.split("=") for el in m.split("; "))
	except:
		logF = open(bn+".log", "a")
		logF.write(line)
		logF.close()
		continue
	l = []
	for w in which:
		if w=="path":
			l.append(f)
			continue
		el = d.get(w, "NA").strip('"')
		el = el if el != "" else "NA"
		l.append(el)
	print "\t".join(l)
input.close()
exit()
