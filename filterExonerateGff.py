#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='Parses gff lines from exonerate output and filter the alignments.')
parser.add_argument('-i', '--input', type=str, default='stdin', help='Input file name. %(default)s')
parser.add_argument('-m', '--method', type=str, help='Method for filtering: <best|all>')
#parser.add_argument('-u', '--up', type=int, default='0', help='How many nucleotides to extend upstream')
#parser.add_argument('-d', '--down', type=int, default='0', help='How many nucleotides to extend downstream')
#parser.add_argument('--starts', action='store_true', default=False, help='Take only the first nucleotide')
args = parser.parse_args()


method = args.method
attr_list = ['gene_id', 'transcript_id']

input = sys.stdin if args.input=="stdin" else open(args.input)

d = {}
gff = False
for line in input:
	if line.startswith("# --- START OF GFF DUMP ---"):
		gff = True
	if line.startswith("# --- END OF GFF DUMP ---"):
		gff = False
	if line.startswith("#") and gff:
		continue
	if gff and "Target" not in line:
		line_sp = line.strip("\n").split("\t")
		chr, ann, el, start, end, score, strand, score2, attr = line_sp
		if el == "similarity":
			continue
		chr, subseq = chr.split(":")
		genomeStart, genomeEnd = map(int, subseq.split("-"))
		#genomeStart, length = map(int, subseq.split("(")[1].strip(")").split(","))
		liftStart = str(int(start)+genomeStart)
		liftEnd = str(int(end)+genomeStart)
		attr = attr.replace(" ;", ";").strip() + ";"
		attrD = dict(a.split(" ") for a in attr.strip(";").split("; ") if attr != ";")
		if el == "gene":
			k = "_".join((chr, start, end, strand, attrD["sequence"]))
			d.setdefault(k, {})
			attrD["gene_id"] = k
			attrD["transcript_id"] = k
			attr = " ".join((a+" \""+attrD[a]+"\";" for a in attr_list))
			attr = attr + " " + " ".join((a+" \""+attrD[a]+"\";" for a in attrD.iterkeys() if a not in attr_list))
			liftLine = "\t".join((chr, ann, el, liftStart, liftEnd, score, strand, score2, attr))
			d[k]["score"] = int(score)
			d[k]["gene"] = liftLine
		if el == "exon":
			score = str(d[k]["score"])
			attrD["gene_id"] = k
			attrD["transcript_id"] = k
			attr = " ".join((a+" \""+attrD[a]+"\";" for a in attr_list))
			attr = attr + " " + " ".join((a+" \""+attrD[a]+"\";" for a in attrD.iterkeys() if a not in attr_list))
			liftLine = "\t".join((chr, ann, el, liftStart, liftEnd, score, strand, score2, attr))
			d[k].setdefault("exons", []).append(liftLine)

if not d:
	print >> sys.stderr, "Empty alignment" 
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


