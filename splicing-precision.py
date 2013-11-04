#!/usr/bin/env python

# Import modules
from argparse import ArgumentParser
from subprocess import Popen, PIPE
import pybedtools as BT


# Argument parsing

parser = ArgumentParser(description='Compute the ratio of split maps spanning annotated junctions vs unannotated nearby junctions.')
parser.add_argument("-s", "--ssj", type=str, help="file with junction read counts.")
parser.add_argument("-g", "--gtf", type=str, help="gtf file with all elements (genes, transcripts, exons).")
parser.add_argument("-w", "--window", type=int, help="number of nucleotides around junctions.", default=20)
args = parser.parse_args()


# BEGIN

gstart, gend = {}, {}
tstart, tend = {}, {}
gexon, texon = {}, {}
gjunction, tjunction = {}, {}

# ==!!!== gtf has offset=1 ==!!!==

for i,line in enumerate(open(args.gtf, "r")):
#	if i > 15000:
#		break
	chr, ann, el, start, end, score1, strand, score2, tag = line.split("\t")
	d = dict(item.split() for item in tag.strip().strip(";").split("; "))
	gene, transcript = d["gene_id"].strip("\""), d["transcript_id"].strip("\"")
	start, end = int(start), int(end)
	exon = "_".join((chr, str(start-1), str(end), strand))
	junctionL = "_".join((chr, str(start-1), str(start), strand))
	junctionR = "_".join((chr, str(end-1), str(end), strand))
	if el == "exon":
		gexon.setdefault(exon, []).append(gene)
		texon.setdefault(exon, []).append(transcript)
		gjunction.setdefault(junctionL, []).append(gene)
		gjunction.setdefault(junctionR, []).append(gene)
		tjunction.setdefault(junctionL, []).append(transcript)
		tjunction.setdefault(junctionR, []).append(transcript)
	if el == "gene":
#		if strand == "-": start, end = end, start
		gstart[gene] = int(start-1)
		gend[gene] = int(end)
	if el == "transcript":
#		if strand == "-": start, end = end, start
		tstart[transcript] = int(start-1)
		tend[transcript] = int(end)
		


def checkJunction(junction):
	chr, start, end, strand = junction.split("_")
	start, end = int(start), int(end)
	genes, transcripts = ",".join(gjunction[junction]), ",".join(tjunction[junction])
	if start not in set(tstart[tx] for tx in tjunction[junction]) and end not in set(tend[tx] for tx in tjunction[junction]):
		return "\t".join((chr, str(end-1), str(end), genes, transcripts, strand))
	return "\n"



ssj = BT.BedTool(args.ssj)
a = BT.BedTool("\n".join(checkJunction(junction) for junction in gjunction.keys()), from_string=True)

a.window(ssj, w=args.window).saveas("test.windowBed")



exit()

