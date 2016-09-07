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

parser.add_argument('-f', '--file', type=str, default="stdin" ,
	help='tsv file with the path of the bigWig (absolute or local) and metadata. [default=stdin]')

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

parser.add_argument("-b", "--browser", default="ucsc", type=str,
	help="Browser for which the urls are written: <ucsc> <igv> [default=%(default)s]]")

parser.add_argument("--composite", action="store_true", default=False,
	help="If the files are part of a composite track [default=%(default)s]]")

parser.add_argument("--groups", type=str, default=None,
	help="Define groups for composite tracks, comma-separated [default=%(default)s]]")

#parser.add_argument("--overlay", type=str, default=None,
#	help="Define groups for overlay tracks, comma-separated [default=%(default)s]]")


# Read arguments
args = parser.parse_args()

name_fields = args.name_fields
color_factors = args.color_factors
palette = args.palette_file
id = args.id
host = args.host
groups = args.groups.split(",") if args.groups else None


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
open_tsv = sys.stdin if args.file == "stdin" else open(args.file)
mdata = readMetadata(open_tsv)

# ---- Colors ----

# Read color palette
color_palette = [hex_to_rgb(line.strip().split("\t")[0]) for line in open(palette)]
# Read the color levels
color_levels = sorted(set(",".join([d.get(color_factor) for color_factor in color_factors.split(",")]) for f,d in mdata.iteritems()))


print

# Print resource lines for igv browser
if args.browser == "igv":
	print "<Resources>"
	for f,d in sorted(mdata.iteritems()):
		url = host + "/" + f.split("/")[-1]
		print "<Resource path=\"%s\"/>" %(url)
	print "</Resources>"
	print
	print "<Panel height=\"6473\" name=\"Panel1457954443035\" width=\"1581\">"



# PARAMETERS

limitMin = "0"
limitMax = "80"
viewLimits = "%s:%s" %(limitMin, limitMax)
type = "bigWig"
visibility = "2"
autoScale = "off"
maxHeightPixels = "30"

if args.composite:
	print "track Composite"
	print "dragAndDrop subTracks"
	print "allButtonPair on"
	print "compositeTrack on"
	print "shortLabel composite"
	print "longLabel composite"
	print "type %s" %(type)
	print "viewLimits %s" %(viewLimits)
	print "autoScale %s" %(autoScale)
	print "maxHeightPixels %s" %(maxHeightPixels)
	print "visibility %s" %(visibility)
	if args.groups:
		for i,group in enumerate(groups):
			print "subGroup%s %s %s \\" %(i+1, group, group)
			gTags = set(list("%s=%s" %(v,v) for value1 in mdata.itervalues() for k,v in value1.iteritems() if k == group))
			print " \\\n".join(gTags) + ";"
		print "sortOrder %s" %(" ".join(["%s=+" %(group) for group in groups]))
	print ""

for f,d in sorted(mdata.iteritems()):
	basename = f.split("/")[-1]
	url = host + "/" + basename 
	name = '"' + ",".join([d.get(name_field).strip("\"") for name_field in name_fields.split(",")]) + '"'
	color_level = ",".join([d.get(color_factor) for color_factor in color_factors.split(",")])
	color = ",".join(color_palette[color_levels.index(color_level)][:3])

# ------ UCSC ------------
	
	if args.browser == "ucsc":
		
		if args.composite:
			print "track %s" %(basename)
			print "bigDataUrl %s" %(url)
			print "shortLabel %s" %(name)
			print "parent Composite"
			print "color %s" %(color)
			print "type %s" %(type)
			if args.groups:
				print "subGroups " + " ".join(["=".join((group, d[group])) for group in groups])
			print ""
			continue
	
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


# ----- IGV ------

	if args.browser == "igv":
		track_line = {
			"name" : name,
			"altColor" : color,
			"color" : color,
			"type" : type,
			"visibility" : visibility,
			"autoScale" : autoScale,
			"maxHeightPixels" : maxHeightPixels,
			"id" : url,
			"displayMode" : "COLLAPSED",
			"featureVisibilityWindow" : "-1",
			"fontSize" : "12",
			"showReference" : "false",
			"visible" : "true"
		}
		dataRange = {
			"maximum" : limitMax,
			"minimum" : limitMin,
			"baseline" : "0",
			"drawBaseline" : "true",
			"flipAxis" : "false",
			"type" : "LINEAR"
		}
# <Track altColor="0,0,178" autoScale="true" color="175,175,175" colorScale="ContinuousColorScale;0.0;142.0;255,255,255;175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="http://genome.crg.es/~abreschi/ENCODE/prCells/bigWig/ENCLB088PVJ_-strand.bigwig" name="ENCLB088PVJ_-strand.bigwig" showReference="false" snpThreshold="0.2" sortable="true" visible="true">
# <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="173.0" minimum="0.0" type="LINEAR"/>
# </Track>
		print "<Track " + " ".join(["%s=\"%s\"" %(item[0], item[1].strip("\"")) for item in track_line.items()]) + ">"
		print "<DataRange " + " ".join(["%s=\"%s\"" %(item[0], item[1].strip("\"")) for item in dataRange.items()]) + "/>"
		print "</Track>"

if args.browser == "igv":
	print "</Panel>"

print

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
