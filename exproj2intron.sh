#!/usr/bin/env bash

file=$1

# the lines in file should be
# ENSG00000000003.10	chrX_99883667_99884983_m


cat $file |\
sed 's/_/\t/g' |\
sort -k 1,1 -k 3,3n |\

awk '
BEGIN{OFS="\t"}

{
s = ($5 == "p" ? "+":"-");
if (gene[$1] == "") {chr[$1]=$2; strand[$1]=s; gene[$1]=$3","$4; next}
if (gene[$1] != "") {gene[$1] = gene[$1](",")$3(",")$4; next}
}

END{
for (g in gene) {
#	print g, gene[g]
	split(gene[g], a, ",")
	for (i in a) {
		if (i == 1 || i == length(a)) {continue}
		if (i%2 == 0) {print g,chr[g]"_"a[i]"_"a[i+1]"_"strand[g]}
		}
	}
}
' 
