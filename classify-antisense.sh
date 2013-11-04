#!/bin/bash

file1=$1
file2=$2

# NOTE:
# file 1 and 2 must be files with exonic projections
# The output will be a list of gene-pairs in opposite orientation and the classification of their relative position

if [ $# != 2 ]; then
echo ""
echo "USAGE: $0 sense.gtf antisense.gtf"
echo "Note: The first field in the gtf tags must be the gene id"
echo ""
exit 1
fi


#file1=~/Documents/db/human/gencode10/Long/MakeSP/gen10.long.makeSP.exonproj.gff
#file2=gen10.long.exonproj.antisense.gff

# =========
# BEGIN
# =========

intersectBed -a $file1 -b $file2 -wo -S |\

awk -v gtf1=$file1 -v gtf2=$file2 'BEGIN{OFS="\t"; FS="\t"
# read start and end for genes in the first exon proj file
while(getline<gtf1>0){
#gene=$10;split(gene,a,"\"");
split($9, tags, "; "); gene=tags[1]; split(gene,a,"\"");

gnstart[a[2]] = (gnstart[a[2]]=="" ? $4 : min(gnstart[a[2]],$4))
gnend[a[2]] = (gnend[a[2]]=="" ? $5 : max(gnend[a[2]],$5))
strand[a[2]] = $7
	}
# read start and end for genes in the second exon proj file
while(getline<gtf2>0){
#gene=$10;split(gene,a,"\"");
split($9, tags, "; "); gene=tags[1]; split(gene,a,"\"");
gnstart[a[2]] = (gnstart[a[2]]=="" ? $4 : min(gnstart[a[2]],$4))
gnend[a[2]] = (gnend[a[2]]=="" ? $5 : max(gnend[a[2]],$5))
strand[a[2]] = $7
	}
}

function max(x,y){return (x>=y) ? x : y;}
function min(x,y){return (x<=y) ? x : y;}

{
# extract gene names
#split($10,a,"\"");split($20,b,"\"");gn1=a[2];gn2=b[2];
split($9,a,"; "); split($18,b,"; "); 
split(a[1],a1,"\""); split(b[1], b1, "\"");
gn1=a1[2];gn2=b1[2];
pairs[gn1,gn2]
# extract the length of the intersection
ovlap = $NF
# take the min start1
if(start1[gn1,gn2]==""){start1[gn1,gn2]=(start1[gn1,gn2]>$4?$4:start1[gn1,gn2])}
# take the min start2
if(start2[gn1,gn2]==""){start2[gn1,gn2]=(start2[gn1,gn2]>$13?$13:start2[gn1,gn2])}
# take the max end1
if(end1[gn1,gn2]==""){end1[gn1,gn2]=(end1[gn1,gn2]<$5?$5:end1[gn1,gn2])}
# take the max end2
if(end2[gn1,gn2]==""){end2[gn1,gn2]=(end2[gn1,gn2]<$14?$14:end2[gn1,gn2])}
}

END{
for (pair in pairs){
	split(pair,ar, SUBSEP)
	gn1_start=gnstart[ar[1]]
	gn2_start=gnstart[ar[2]]
	gn1_end=gnend[ar[1]]
	gn2_end=gnend[ar[2]]
	strand1=strand[ar[1]]
	strand2=strand[ar[2]]

# debug
#	if (ar[1]=="ENSG00000162385.6") { print "A", gn1_start, gn1_end, gn2_start, gn2_end}


# -- Conditions on relative posisions --

	if (strand1 == "+" && strand2 == "-") { 
# 5 head to head
		if (gn2_start<gn1_start && gn1_start<gn2_end && gn1_end>gn2_end) {class[pair]="5head-to-head"}
# 3 tail to tail
		if (gn2_end>gn1_end && gn2_start>gn1_start && gn2_start<=gn1_end) {class[pair]="3tail-to-tail"}
		}

	if (strand1 == "-" && strand2 == "+") { 
# 5 head to head
		if (gn2_start<gn1_start && gn1_start<=gn2_end && gn1_end>gn2_end) {class[pair]="3tail-to-tail"}
# 3 tail to tail
		if (gn2_end>gn1_end && gn2_start>gn1_start && gn2_start<=gn1_end) {class[pair]="5head-to-head"}
		}


# Internal intronic
	if (gn2_start>=gn1_start && gn2_start<gn1_end && gn2_end>gn1_start && gn2_end<=gn1_end) {class[pair]="internal"}
# External exonic
	if (gn2_start<=gn1_start && gn2_end>=gn1_end) {class[pair]="external"}

	print ar[1],gn1_start,gn1_end,strand1, ar[2],gn2_start,gn2_end,strand2, class[pair]
	}
}'












# *** Schema ***

# Note: there is at least 1 exonic nt overlap


# #####################
# ## 5' head to head ##
# #####################

#                               |---> gn1
#                             5'|===============-------=========--------------========
# 3' ==========-------------========|
#                           gn2 <---|


# #####################
# ## 3' tail to tail ##
# #####################


#     |---> gn1
#   5'|===============-------=========--------------========
#                                                       3' ==========-------------========|
#                                                                                gn2  <---|
#

# #####################
# ## internal        ##
# #####################


#     |---> gn1
#   5'|===============-------=========--------------========
#                 3' ==========-------------========|
#                                          gn2  <---|


# #####################
# ## external        ##
# #####################


#                  |---> gn1
#                5'|===============-------=========--------------========
#  3' =================--------------==========-------------------------------========|
#                                                                            gn2  <---|


