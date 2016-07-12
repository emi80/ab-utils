#!/usr/bin/env python

import sys, argparse, operator


############  ARGPARSE ##############

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='From RAMPAGE bed (after RAMPAGE2gene.py) for different samples, make a matrix')

parser.add_argument('-i', '--input', type=str, default="stdin",
	help='Input file. 2 columns. col1: bed file, col2: sample id (will be the column header) [default=%(default)s]')

parser.add_argument('-M', '--merged', type=str,
	help='File where annotation of merged TSSs is reported (optional)')

# Read arguments
args = parser.parse_args()

def getOverlap(a, b):
	# a and b are tuples of 4 elements: (chr, start, end, strand)
	if len(a) == 0 or len(b) == 0:
		return -1
	return max(0, min(int(a[2]), int(b[2])) - max(int(a[1]), int(b[1])))

def mergeIntervals(a, b):
	if len(a) == 0 and len(b) == 0:
		print "ERROR: Merging empty intervals"
		exit()
	if len(a) == 0:
		return b
	if len(b) == 0:
		return a
	return a[0], str(min(int(a[1]), int(b[1]))), str(max(int(a[2]), int(b[2]))), a[3]

inF = sys.stdin if args.input == "stdin" else open(args.input)


geneDict = dict()
tssDict = dict()

chrs = set()
counter = 1

if args.merged:
	outF = open(args.merged, "w")
	outF.close()

samples = set()
for line in inF:
	bed, id = line.strip().split("\t")
	samples.add(id)
	for line in open(bed):
		tss_chr, tss_start, tss_end, tss_id, score, tss_strand, gene_id, dist, nt, reads, ann = line.strip().split("\t")
		coord = (tss_chr, tss_start, tss_end, tss_strand)
		if ann == "novel":
			gene_id = tss_chr
			chrs.add(tss_chr)
		geneDict.setdefault(gene_id, [[], []])
		geneDict[gene_id][0].append(coord)
		geneDict[gene_id][1].append({id:reads})

print "\t".join(sorted(samples))
for gene_id, (coords, readsDicts) in geneDict.iteritems():
	if len(coords) == 1:
		tssDict[counter] = [coords[0], coords, readsDicts]

		# output
		pos, all_pos, expr = tssDict[counter]
		a = dict(pair for d in expr for pair in d.items())
		line = str(counter) + "_" + "_".join(pos) + "\t" + "\t".join(a.get(id, "0") for id in sorted(samples))
		print line
		if args.merged:
			outF = open(args.merged, "a")
			outF.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(
				str(counter), 
				gene_id,
				"_".join(pos),
				",".join("_".join(el) for el in coords),
				",".join(a.keys()),
				",".join(a.values())
				)
			)
			outF.close()

		counter += 1
		continue
	if len(coords) > 1:
		sortIndexes, sortCoords = zip(*sorted(enumerate(coords), key=lambda x: (x[1][1], x[1][2])))
		merged, expr, tss_merged = (), [], []
		for i, coord in zip(sortIndexes, sortCoords):
			overlap = getOverlap(merged, coord)
			if overlap != 0:
				merged = mergeIntervals(coord, merged) 
				expr.append(readsDicts[i])
				tss_merged.append(coord)
			if overlap == 0:
				tssDict[counter] = [merged, tss_merged, expr]

				# Output
				pos, all_pos, expr = tssDict[counter]
				a = dict(pair for d in expr for pair in d.items())
				line = str(counter) + "_" + "_".join(pos) + "\t" + "\t".join(a.get(id, "0") for id in sorted(samples))
				print line
				if args.merged:
					outF = open(args.merged, "a")
					outF.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(
						str(counter), 
						gene_id,
						"_".join(pos),
						",".join("_".join(el) for el in all_pos),
						",".join(a.keys()),
						",".join(a.values())
						)
					)
					outF.close()

				counter += 1
				merged, expr, tss_merged = coord, [readsDicts[i]], [coord]
		tssDict[counter] = [merged, tss_merged, expr]

		# Output
		pos, all_pos, expr = tssDict[counter]
		a = dict(pair for d in expr for pair in d.items())
		line = str(counter) + "_" + "_".join(pos) + "\t" + "\t".join(a.get(id, "0") for id in sorted(samples))
		print line
		if args.merged:
			outF = open(args.merged, "a")
			outF.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(
				str(counter), 
				gene_id,
				"_".join(pos),
				",".join("_".join(el) for el in all_pos),
				",".join(a.keys()),
				",".join(a.values())
				)
			)
			outF.close()

		counter += 1
		merged, expr, tss_merged = (), [], []
		#print tssDict
		#break


#print "\t".join(sorted(samples))
#for tss, (pos, coords, expr) in tssDict.iteritems():
#	a = dict(pair for d in expr for pair in d.items())
#	line = str(tss) + "_" + "_".join(pos) + "\t" + "\t".join(a.get(id, "0") for id in sorted(samples))
#	print line

exit()
