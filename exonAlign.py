#!/usr/bin/env python

import argparse, sys
import subprocess as sp
import itertools as it


parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Align......')
parser.add_argument("-a", "--file1", type=str, help="file name or stdin")
parser.add_argument("-b", "--file2", type=str, help="file name or stdin")
#parser.add_argument("-x", type=str, default="1", help="index of the field you want file a to be joined on. Accept multiple indeces comma-separated")
#parser.add_argument("-y", type=str, default="1", help="index of the field you want file b to be joined on. Accept multiple indeces comma-separated")
#parser.add_argument("--a_header", action="store_true", help="file \"a\" has a header")
#parser.add_argument("--b_header", action="store_true", help="file \"b\" has a header")
parser.add_argument("--ali", help="alignment file, fasta format")
args = parser.parse_args()


def read_gtf(f, outF):
	d = {}
	for line in f:
		chr, ann, element, start, end, score1, strand, score2, tag = line.strip().strip(";").split("\t")
		if element != "CDS":
			continue
		tag_d = dict(el.split(" ") for el in tag.split("; "))
		gene_id = tag_d["transcript_id"].strip("\"")
		d.setdefault(gene_id,{})
		i = max(d[gene_id].keys() + [0]) + 1
		d[gene_id][i] = (start, end, strand)
	print d
	for gene_id in d.iterkeys():
		lastExon = max(d[gene_id].keys())
		gene_end = int(d[gene_id][lastExon][1])
		gene_start = int(d[gene_id][1][0])
		for i in sorted(d[gene_id].keys()):
			# Write the exons
			start = int(d[gene_id][i][0])
			end = int(d[gene_id][i][1])
			if strand == "+":
				start = start - gene_start
				end = end - gene_start
			if strand == "-":
				start = gene_end - start# + gene_start
				end = gene_end - end# + gene_start
			lineOut = "\t".join((gene_id, str(start), str(end), "exon%s" %i, "struct"))
			outF.write(lineOut + "\n")
			print lineOut
	
			# Write the introns	
			if i == lastExon:
				continue
			start = int(d[gene_id][i][1]) + 1
			end = int(d[gene_id][i+1][0]) - 1
			if strand == "+":
				start = start - gene_start
				end = end - gene_start
			if strand == "-":
				start = gene_end - start# + gene_start 
				end = gene_end - end# + gene_start
			lineOut = "\t".join((gene_id, str(start), str(end), "intron%s" %i, "struct"))
			outF.write(lineOut + "\n")
			print lineOut

	return d


out = "exonAlign.out.tsv"
#outF = open(out, "w")
outF = open(out, "w")

gtf = read_gtf(sp.Popen("sort -k 1,1 -k 4,4n " + args.file1, shell=True, stdout=sp.PIPE).stdout, outF)
if args.file2:
	gtf.update(read_gtf(sp.Popen("sort -k 4,4n " + args.file2, shell=True, stdout=sp.PIPE).stdout, outF))
outF.close()


cod = 3

d = {}
unionBound = set()
indexBound = {}
names = list()

aln = open(args.ali).readlines()
for NR,line in enumerate(aln):
	line = line.strip()
	if line == "":
		continue
	if line.startswith(">"):
		seq = ""
		name = line.strip(">").split()[0]
		names.append(name)
		d[name] = {}
		strand = gtf[name][1][2]
		# Set with the all exon boundaries
		boundaries = set(int(v[i]) for v in gtf[name].values() for i in (0,1))
		indexBound[name] = set()
		if strand == "+":
			firstExon = 1
			lastExon = max(gtf[name].keys())
			gene_start = int(gtf[name][firstExon][0])
			gene_end = int(gtf[name][lastExon][1])
		if strand == "-":
			firstExon = max(gtf[name].keys())
			lastExon = 1
#			gene_start = int(gtf[name][lastExon][1])
#			gene_end = int(gtf[name][firstExon][0])
			gene_start = int(gtf[name][lastExon][0])
			gene_end = int(gtf[name][firstExon][1])
		continue
	seq = seq + line
	if NR+1 != len(aln):
		if not aln[NR+1].startswith(">"):
			continue
	for i,l in enumerate(list(seq)):
		if l == "-":
			continue
		if not d[name]:
			exon = firstExon
			offset = 0
		if strand == "+":
			for slide in range(cod):
				pos = int(gtf[name][exon][0]) + offset
				d[name].setdefault(i,[]).append((pos - gene_start, exon))
#				print name, d[name]
				if pos in boundaries: indexBound[name].add(i)
				print name,i,l,pos,exon
				offset = offset + 1 
				if pos + 1 > int(gtf[name][exon][1]):
					exon = exon + 1
					offset = 0
		if strand == "-":
			for slide in range(cod):
				pos = int(gtf[name][exon][1]) - offset
				flip_pos = gene_end - pos #+ gene_start
				d[name].setdefault(i,[]).append((flip_pos, exon))
				### if in boundaries: add
				if pos in boundaries: indexBound[name].add(i)
				print name,i,l,pos,exon
#				print name, d[name]
				offset = offset + 1 
				if pos - 1 < int(gtf[name][exon][0]):
					print pos, offset, gtf[name]
					exon = exon - 1
					offset = 0

unionBound = set(el for s in indexBound.values() for el in s)
print unionBound

outF = open(out, "a")
#for gene1, gene2 in it.combinations(names, 2):
for i in range(len(names)-1):
	gene1,gene2 = names[i:i+2]
	# Get all the boundaries between the two genes
	pairBound = indexBound[gene1].union(indexBound[gene2])
	intersectBound = unionBound.intersection(pairBound)
	for i in intersectBound:
		for slide in range(cod):
			try:
				edge1 = d[gene1][i][slide][0]
			except KeyError:
				edge1 = "NA"
			try: 
				edge2 = d[gene2][i][slide][0]
			except KeyError:
				edge2 = "NA"
			lineOut = "\t".join((gene1, str(edge1), str(edge2), gene2, "conn"))
			outF.write(lineOut + "\n")
			print lineOut
#	for bound in intersectBound:
#		print

