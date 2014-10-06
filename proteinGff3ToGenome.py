#!/usr/bin/env python

import sys
import argparse
sys.path.append('/users/rg/mmariotti/Scripts')
from MMlib import gene

############  ARGPARSE ##############

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='Print track lines for a set of bigWig from a tsv file')

parser.add_argument('-f', '--file', type=str, 
	help='gff3 with protein coordinates. Can be stdin')

parser.add_argument("-b", "--bed12", type=str, 
	help="bed12 with the CDS")

parser.add_argument("-o", "--output", type=str, default="stdout",
	help="name of the output file. [default=%(default)s]")

parser.add_argument("-u", "--ucsc", type=str,
	help="Specify a name for an optional file to be loaded as UCSC tracks")

parser.add_argument('-n', "--name_fields", type=str, 
	help="gff3 keys, whose values are pasted to give the sequence name (comma-separated).")

#parser.add_argument("-i", "--id", type=str, 
#	help="name of the fields which contains the path (can be local or full)")
#
#parser.add_argument("-H", "--host", type=str, 
#	help="internet directory that contains the tracks")
#
#parser.add_argument("-c", "--color_factors", default="cell", type=str, 
#	help="comma-separated fields that are pasted to be the color factor. [default=%(default)s]")
#
#parser.add_argument("-p", "--palette_file", default="/users/rg/abreschi/R/palettes/Accent.8.txt", type=str, 
#	help="File with the color palette in hex format. Transparency is ignored. [default=%(default)s]")

# Read arguments
args = parser.parse_args()

#name_fields = args.name_fields
#color_factors = args.color_factors
#palette = args.palette_file
#id = args.id
#host = args.host


# FUNCTIONS

# Read the bed12
def readBed12(bed12):
	ref = {}
	# NB: Marco's lib uses coordinates with offset 1
	for line in open(bed12, 'r'):
#		chr, start, end, id, score, strand, Tstart, Tend, itemRgb, blockCount, blockSizes, blockStarts 
		line_el = line.strip().split("\t")
		chr, start, end, id, score, strand = line_el[0:6]
		Tstart, Tend = line_el[6:8]
		start, end, Tstart, Tend = map(int, (start, end, Tstart, Tend))
		blockCount, blockSizes, blockStarts = line_el[9:12]
		blockSizes = map(int, blockSizes.strip(",").split(","))
		blockStarts = map(int, blockStarts.strip(",").split(","))
		g = gene()
		g.chr = chr
		g.strand = strand
		for i in range(int(blockCount)):
			if i == 0:
				CDS = [Tstart, start+blockStarts[i]+blockSizes[i]]
			else:
				CDS = [start+blockStarts[i], start+blockStarts[i]+blockSizes[i]]
			g.add_exon(CDS[0], CDS[1])
		ref[id] = g
	return ref


# Convert hexadecimal colors to rgb
def hex_to_rgb(value):
    value = value.lstrip('#').lower()
    lv = len(value)
    return tuple(str(int(value[i:i+lv/3], 16)) for i in range(0, lv, lv/3))		


#~~~~~~~~~~#
# BEGIN    #
#~~~~~~~~~~#

# Read input
input = sys.stdin if args.file == 'stdin' else open(args.file, 'r')

# Read output file
output = sys.stdout if args.output == "stdout" else open(args.output, "w")

parents = readBed12(args.bed12)

# Read the keys to use for description
fields = args.name_fields

# Palette
palette = "/users/rg/abreschi/R/palettes/Paired.12.txt"
color_palette = [",".join(hex_to_rgb(line.strip())) for line in open(palette)]


# Read gff3 and convert the coordinates
tracks = {}
for line in input:
	if line.startswith("##"):
		continue
	if len(line.split("\t")) != 9:
		continue
	chr, ann, element, start, end, score1, strand, score2, tags = line.strip().strip(";").split("\t")
	if ann == ".":
		continue
	tags_d = dict(item.split("=") for item in tags.split(";"))
	desc = "_".join(tags_d.get(k, "NA") for k in fields.split(","))
	g = gene()
	g.chr = chr
	g.add_exon(int(start) * 3 +1, int(end) * 3 +1)
	g.strand = strand
	parent = parents[chr]
	g.restore_absolute_coordinates(parent)
	attrs = 'gene_id "%s"; transcript_id "%s";' %(chr, desc)
	for exon in g.exons:
		out_line = "\t".join((
			parent.chr, 
			ann,
			"exon", 
			str(exon[0]),
			str(exon[1]), 
			".", 
			parent.strand, 
			".", 
			attrs
		))
		output.write(out_line+"\n")
	if args.ucsc is not None:
		tracks.setdefault(ann, []).append(out_line+"\n")

# TODO: different colors for different tracks :)

output.close()
	
# Write the UCSC track file if requested		
if args.ucsc is not None:
	output_tracks = open(args.ucsc, "w")
	for i,key in enumerate(tracks.iterkeys()):
		output_tracks.write("track name=%s type=gtf color=%s\n" %(key, color_palette[i]))
		for value in tracks.get(key):
			output_tracks.write(value)
	output_tracks.close()

exit()	


