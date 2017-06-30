#!/usr/bin/env python

import sys, argparse, gzip
import subprocess as sp
from operator import itemgetter


############  ARGPARSE ##############

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='Assign a gene id to RAMPAGE output')

parser.add_argument('-R', '--RAMPAGEgff', type=str, default="stdin",
	help='RAMPAGE output, gff [default=%(default)s]')

#parser.add_argument('-o', '--output', type=str, default="stdout",
#	help='Output file name. [default=%(default)s]')

parser.add_argument('-a', '--annotation', type=str,
	help='Annotation file (what matters here is transcript coordinates)')


# Read arguments
args = parser.parse_args()

inF = open(args.RAMPAGEgff)
if args.RAMPAGEgff.endswith(".gz"):
	inF = gzip.open(args.RAMPAGEgff)
if args.RAMPAGEgff == "stdin":
	inF = sys.stdin


tssDict = dict()
reads = dict()

p = sp.Popen('intersectBed -a - -b %s -s -wao' %(args.annotation), shell=True, stdin=sp.PIPE, stdout=sp.PIPE)

bedLines = ""
for line in inF:
	if line.startswith("track"):
		continue
	if line.strip() == "":
		continue
	chr, ann, el, start, end, score1, strand, score2, tags = line.strip().split("\t")
	tags = dict(el.split(" ") for el in tags.strip(";").split("; "))
	g_chr, g_strand, g_start, g_end = tags["gene_id"].strip("\"'").split("_")
	
	tss_id = tags["tss_id"].strip("\"'")

	# for each tss id store its coordinates
	tssDict[tss_id] = (chr, start, end, strand)

	reads[tss_id] = score1

	# Bed file with the coordinates of transcribed regions (not the tss!)
	bedLine = "\t".join((g_chr, g_start, g_end, tags["tss_id"].strip("\"'"), "0", strand))
#	p.stdin.write(bedLine + "\n")
	bedLines = bedLines + bedLine + "\n"
#p.stdin.close()
#print "done"

overlap, distance = dict(), dict()
novel = set()

#chr11	GRIT	TSS	62609279	62609284	766	-	.	gene_id 'chr11_minus_62599813_63183458'; gene_name 'chr11_minus_62599813_63183458'; tss_id 'TSS_chr11_minus_62599813_63183458_pk1'; peak_cov '5,11,0,0,750';


for line in p.communicate(input = bedLines)[0].strip().split("\n"):
	# Parse the output of intersectBed
	bed, gff, nt = line.split("\t")[:6], line.strip().split("\t")[6:-1], int(line.strip().split("\t")[-1])
	tss_id = bed[3]
	tss_chr, tss_start, tss_end, tss_strand = tssDict[tss_id]

	if nt != 0:
		tags = dict(el.split(" ") for el in gff[8].strip(";").split("; "))
		gene_id = tags["gene_id"].strip("\"'")
		g_chr, g_ann, g_el, g_start, g_end, g_score1, g_strand = gff[:7]

		# Overlap of transcribed region with annotated gene id
		overlap.setdefault(tss_id, {}).setdefault(gene_id, nt)
		overlap[tss_id][gene_id] = max(overlap[tss_id][gene_id], nt)

		# Distance of TSS to the beginning of the overlapping gene
		if tss_strand == "+":
			dist = abs(int(tss_start) - int(g_start))
		if tss_strand == "-":
			dist = abs(int(tss_end) - int(g_end))
		distance.setdefault(tss_id, {}).setdefault(gene_id, dist)
		distance[tss_id][gene_id] = min(distance[tss_id][gene_id], dist)

	# No gene overlapping
	if nt == 0:
		novel.add(tss_id)


for tss_id in set(tssDict.keys()) - novel:

	# Find closest TSSs (in case of ties report all)
	closest_v = min(distance[tss_id].values())
	closest_k = [k for k in distance[tss_id] if distance[tss_id][k] == closest_v]
	dist = closest_v

	# If the minimum distance to annotated TSSs is larger than 10000, assign the TSS to novel	
	if dist > 10000:
		novel.add(tss_id)
		continue

	if len(closest_k) == 1:
		gene_id = closest_k[0]

	# If more genes have the same distance, take the one with the highest overlap
	if len(closest_k) > 1:
		d = dict((k, overlap[tss_id][k]) for k in closest_k)
		longest_v = max(d.values())
		longest_k = [k for k in d if d[k] == longest_v]
		# If more genes have the same exact overlap take the first alphabetically
		gene_id = sorted(longest_k)[0]

	nt = overlap[tss_id][gene_id]
	tss_chr, tss_start, tss_end, tss_strand = tssDict[tss_id]
	line = "\t".join((tss_chr, tss_start, tss_end, tss_id, str(0), tss_strand, gene_id, str(dist), str(nt), reads[tss_id], "linked"))
	print line

for tss_id in novel:
	tss_chr, tss_start, tss_end, tss_strand = tssDict[tss_id]
	line = "\t".join((tss_chr, tss_start, tss_end, tss_id, str(0), tss_strand, "NA", str(-1), str(-1), reads[tss_id], "novel"))
	print line
	

exit()	
#for line in p.stdout:
#	print line.strip()
#	break
