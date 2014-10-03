#!/usr/bin/env python

import sys
import argparse

#host = "http://genome.crg.es/~abreschi/ERC/dmel/UCSC_tracks/"
#
#id = "path"
#
#name_fields = "cell,tissue,view"
#color_factors = "view"
#palette = "/users/rg/abreschi/R/palettes/rainbow.5.txt"

############  ARGPARSE ##############

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='Print track lines for a set of bigWig from a tsv file')

parser.add_argument('-f', '--file', type=str, 
	help='tsv file with the path of the bigWig (absolute or local) and metadata. Can be stdin')

parser.add_argument("-i", "--id", type=str, 
	help="name of the fields which contains the path (can be local or full)")

parser.add_argument("-H", "--host", type=str, 
	help="internet directory that contains the tracks")

parser.add_argument('-n', "--name_fields", default="cell", type=str, 
	help="comma-separated fields that are pasted to give the track name. [default=%(default)s]")

parser.add_argument("-c", "--color_factors", default="cell", type=str, 
	help="comma-separated fields that are pasted to be the color factor. [default=%(default)s]")

parser.add_argument("-p", "--palette_file", default="/users/rg/abreschi/R/palettes/Accent.8.txt", type=str, 
	help="File with the color palette in hex format. Transparency is ignored. [default=%(default)s]")

# Read arguments
args = parser.parse_args()

name_fields = args.name_fields
color_factors = args.color_factors
palette = args.palette_file
id = args.id
host = args.host

##### FUNCTIONS ############

def readIndex(index):
	mdata = {}
	for line in index:
		path, data = line.strip().split("\t")
		d = dict(kv.split("=") for kv in data.strip(";").split("; "))
		mdata[d.get(id)] = d
	return mdata


def readMetadata(tsv):
	mdata = {}
	for i,line in enumerate(tsv):
		if i == 0:
			header = line.strip().split("\t")
			continue
		values = line.strip().split("\t")
		mdata[values[header.index(id)]] = dict((header[i],v) for i,v in enumerate(values))
	return mdata

def hex_to_rgb(value):
    value = value.lstrip('#').lower()
    lv = len(value)
    return tuple(str(int(value[i:i+lv/3], 16)) for i in range(0, lv, lv/3))		



#mdata = readMetadata("/users/rg/abreschi/Documents/ERC/dmel/erc_dmel_rna_dashboard.tsv")


# Read the tsv
open_tsv = sys.stdin if args.file == "stdin" else  open(args.file)
mdata = readMetadata(open_tsv)

# ---- Colors ----

# Read color palette
color_palette = [hex_to_rgb(line.strip()) for line in open(palette)]
# Read the color levels
color_levels = sorted(set(",".join([d.get(color_factor) for color_factor in color_factors.split(",")]) for f,d in mdata.iteritems()))


for f,d in mdata.iteritems():
	url = host + f.split("/")[-1]
	name = '"'+ ",".join([d.get(name_field).strip("\"") for name_field in name_fields.split(",")]) + '"'
	color_level = ",".join([d.get(color_factor) for color_factor in color_factors.split(",")])
	color = ",".join(color_palette[color_levels.index(color_level)][:3])
	viewLimits = "0:80"
	type = "bigWig"
	visibility = "2"
	autoScale = "off"
	maxHeightPixels = "35"
	track_line = {
		"name" : name,
		"viewLimits" : viewLimits,
		"color" : color,
		"type" : type,
		"visibility" : visibility,
		"autoScale" : autoScale,
		"maxHeightPixels" : maxHeightPixels,
		"bigDataUrl" : url
	}
	print "track " + " ".join(["=".join(item) for item in track_line.items()])

exit()

#ls *bw | awk  'BEGIN{print}
#{viewLimits="0:30"; color="0,255,0"
#if ($0~/DNA_methyl/) {viewLimits="0:1"; color="255,127,36"}
#if ($0~/RNAseq/) {viewLimits="0:50"; color="100,149,237"}
#if ($0~/H3K27me3/) {viewLimits="0:5"; color="0,0,0"}
#if ($0~/H3K36/) {viewLimits="0:10"; color="69,139,0"}
#if ($0~/H3K27ac/) {viewLimits="0:10"; color="34,139,34"}
#if ($0~/H3K9me3/) {color="139,0,139"}
#if ($0~/Input/) {color="169,169,169"}
#print "track name="$0" type=bigWig visibility=2 autoScale=off maxHeightPixels=30 color="color" viewLimits="viewLimits" bigDataUrl=http://genome.crg.es/~abreschi/BluePrint/UCSC_tracks/"$0}' 
#
