#!/usr/bin/env python

import random
random.seed(123)
import sys
import argparse
sys.path.append('/users/rg/mmariotti/Scripts')
from MMlib import gene

############  ARGPARSE ##############

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-f', '--file', type=str, default="stdin", 
	help='gff3 with protein coordinates [default=stdin]')

parser.add_argument("-b", "--bed12", type=str, 
	help="bed12 with the CDS. 4th col should be the chr of the input file")

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
		# No CDS in the transcript
		if Tstart == Tend:
			continue
		start, end, Tstart, Tend = map(int, (start, end, Tstart, Tend))
		blockCount, blockSizes, blockStarts = line_el[9:12]
		blockSizes = map(int, blockSizes.strip(",").split(","))
		blockStarts = map(int, blockStarts.strip(",").split(","))
		g = gene()
		g.chromosome = chr
		g.strand = strand
		for i in range(int(blockCount)):
			blockStart = start + blockStarts[i]
			blockEnd = start + blockStarts[i] + blockSizes[i]
			if blockStart <= Tend and blockEnd >= Tstart:
				CDS = [max(Tstart, blockStart), min(Tend, blockEnd)]
				g.add_exon(CDS[0]+1, CDS[1])
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

# Dictionary with tanscript ids and MMlib gene object
parents = readBed12(args.bed12)

# Read the keys to use for description
fields = args.name_fields

# Palette
palette = "/users/rg/abreschi/R/palettes/Paired.12.txt"
color_palette = [",".join(hex_to_rgb(line.strip())) for line in open(palette)]

convert = True 
# Read gff3 and convert the coordinates
if convert:
	tracks = {}
	for line in input:
		if line.startswith("##"):
			continue
		if len(line.split("\t")) != 9:
			continue
		chr, ann, element, start, end, score1, strand, score2, tags = line.strip().strip(";").split("\t")
		if ann == ".":
			continue
		if "?" in start or "?" in end:              # Sometimes there can be a question mark in uniprot position
			continue
		if "<" in start or ">" in end:              # Sometimes there can be a uncertain sign mark in uniprot position
			continue
		tags_d = dict(item.split("=") for item in tags.split(";"))
		desc = "_".join(tags_d.get(k, "NA") for k in fields.split(","))
		g = gene()
		g.chromosome = chr
		g.add_exon(int(start) * 3 -2, int(end) * 3)
		g.strand = parents[chr].strand
		parent = parents[chr]
		try:
			g.restore_absolute_coordinates(parent)
		except:
			print >> sys.stderr, 'WARNING: %s length shorter than protein sequence' %chr
		# chr: the sequence id where the domain hit is found (should correspond to 4th col in bed12 
		# desc: the concatenation of several fields in the gff input
		attrs = 'gene_id "%s"; transcript_id "%s_%s";' %(chr, chr, desc)
#		attrs = 'gene_id "%s"; transcript_id "%s";' %(chr, desc)
		for exon in g.exons:
			out_line = "\t".join((
				parent.chromosome, 
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
		try:
			color_i = color_palette[i]
		except IndexError:
			color_i = ",".join(map(lambda x: str(random.randint(0, 255)), range(3)))
		output_tracks.write("track name=%s type=gtf color=%s\n" %(key, color_i))
		for value in tracks.get(key):
			output_tracks.write(value)
	output_tracks.close()

exit()	


