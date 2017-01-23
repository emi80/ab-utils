#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='Parses gff lines from exonerate output and filter the alignments.')
parser.add_argument('-i', '--input', type=str, default='stdin', help='Input file name. %(default)s')
parser.add_argument('-m', '--method', type=str, help='Method for filtering: <best|all>')
parser.add_argument("--query", action="store_true", default=False, help="Write only the query gff")
#parser.add_argument('-u', '--up', type=int, default='0', help='How many nucleotides to extend upstream')
#parser.add_argument('-d', '--down', type=int, default='0', help='How many nucleotides to extend downstream')
#parser.add_argument('--starts', action='store_true', default=False, help='Take only the first nucleotide')
args = parser.parse_args()


method = args.method
attr_list = ['gene_id', 'transcript_id']

input = sys.stdin if args.input=="stdin" else open(args.input)

d = {}
t = {}
gff, query_gff, target_gff = False, False, False
target_counter, query_counter = 0, 0

for line in input:
	if line.startswith("# --- START OF GFF DUMP ---"):
		gff = True
	if line.startswith("# --- END OF GFF DUMP ---"):
		gff, query_gff, target_gff = False, False, False
	if line.startswith("#") and gff:
		continue
	if gff:
		line_sp = line.strip("\n").split("\t")
		chr, ann, el, start, end, score, strand, score2, attr = line_sp
		attr = attr.replace(" ;", ";").strip() + ";"
	if gff and "Target" in line:
		query_gff = True
		query_counter = query_counter + 1
	if gff and "sequence" in line:
		target_gff = True
		target_counter = target_counter + 1 

# Parse the gff lines of target
	if gff and target_gff:
		if el == "similarity":
			continue
	#	print line
		if ":" in chr:
			chr, subseq = chr.split(":")
			genomeStart, genomeEnd = map(int, subseq.split("-"))
		else:
			genomeStart = 0
		liftStart = str(int(start)+genomeStart)
		liftEnd = str(int(end)+genomeStart)
		attrD = dict(a.split(" ") for a in attr.strip(";").split("; ") if attr != ";")
		if el == "gene":
			#k = "_".join((chr, start, end, strand, attrD["sequence"]))
			k = "%s_match%s" %(attrD["sequence"], target_counter)
			d.setdefault(k, {})
			attrD["gene_id"] = k
			attrD["transcript_id"] = k
			attrD["score"] = score
			attr = " ".join((a+" \""+attrD[a]+"\";" for a in attr_list))
			attr = attr + " " + " ".join((a+" \""+attrD[a]+"\";" for a in attrD.iterkeys() if a not in attr_list))
			liftLine = "\t".join((chr, ann, el, liftStart, liftEnd, ".", strand, score2, attr))
			d[k]["score"] = int(score)
			d[k]["gene"] = liftLine
		if el == "exon":
			#score = str(d[k]["score"])
			attrD["gene_id"] = k
			attrD["transcript_id"] = k
			attrD["score"] = str(d[k]["score"])
			attr = " ".join((a+" \""+attrD[a]+"\";" for a in attr_list))
			attr = attr + " " + " ".join((a+" \""+attrD[a]+"\";" for a in attrD.iterkeys() if a not in attr_list))
			liftLine = "\t".join((chr, ann, el, liftStart, liftEnd, ".", strand, score2, attr))
			d[k].setdefault("exons", []).append(liftLine)


# Parse the gff lines of query
	if gff and query_gff:
		attrD = {}
		attrD["ID"] = "%s_match%s" %(chr, query_counter)
		attrD["score"] = score
		for a in attr.strip(";").split("; "):
			k,v = a.split(" ")[0], a.split(" ")[1:]
			if k == "Align":
				attrD.setdefault(k, []).append(v)
			else:
				attrD[k] = ",".join(v)
		for align in attrD["Align"]:
			qStart, tStart, length = map(int, align)
			qAttr = ";".join(("%s=%s" %(k,v) for k,v in attrD.iteritems() if k != "Align"))
			qLine = "\t".join((chr, ann, el, str(qStart), str(qStart+length), ".", strand, score2, qAttr))
			if args.query:
				print qLine

if args.query:				
	exit()


if not d:
	print >> sys.stderr, "Empty alignment: %s" %args.input
	exit()

if method == "best":
	k,v = max(d.iteritems(), key=lambda (k,v):(v["score"]))
	print v["gene"].strip()
	print "\n".join(v["exons"])

if method == "all":
	for v in d.itervalues():
		print v["gene"].strip()
		print "\n".join(v["exons"])


exit()


