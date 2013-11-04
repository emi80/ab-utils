#!/usr/bin/env bash

# This script adds gtf tags from a tab-separated file based on the gene name
# *** Sep. 5th, 2013 ***

if [ $# != 3 ]; then
echo "USAGE: $0 file.gtf file.tsv index"
echo 
echo "Note: accepts also gene name without '.'"
exit 1
fi

# Arguments

gtf=$1               # gtf file you want to add the tags to (gene gtf)
tsv=$2               # file with values for tags, header is the key, the first column is gene name
fields=$3            # indexes of the columns you want to include in the gtf


# BEGIN

awk -v tsv=$tsv -v fields=$fields '
BEGIN{
OFS="\t";
split(fields, selFields, ",");
i=1
while(getline<tsv>0) {
	if (i==1) {for (f in selFields) {header[selFields[f]]=$selFields[f]}; i++; continue}
	for (f in selFields) {toAdd[$1, header[selFields[f]]] = $selFields[f]}
	}
}

{
gene=$10; split(gene, a,"\""); split(a[2], b, ".");
printf $0;
for (f in selFields) {
	key = header[selFields[f]]
	# if the value is not found add write NA. (or an empty string is better?)
	value = ( toAdd[a[2], key] == "" ? (toAdd[b[1], key] == "" ? "NA" : toAdd[b[1], key]) : toAdd[a[2], key] )
	printf " "key" \""value"\";"
	}
printf "\n"
}

END{}
' $gtf











# Examples of input

# -- gtf
# chr1	HAVANA	gene	11869	14412	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";

# -- tsv
# gene	EIR	JR	nNSP	nSPL	coSI
# ENSG00000254659.1	0	0	6	6	NA
# ENSG00000132323.4	71	226	18	26	0.815067

# -- fields
# 2,3,4
