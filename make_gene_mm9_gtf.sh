#!/bin/bash

cat - |\

awk '
{
genes[$10] = $10
if ($1 !~/^Contig/ && $1 !~/^Ultra/) {
	chr[$10] = "chr"$1
} else {chr[$10] = $1}
if (start[$10] == "" || start[$10] > $4 ) {start[$10] = $4}
if (end[$10] == "" || end[$10] < $5) {end[$10] = $5}
strand[$10] = $7
gene_type[$10] = $18
gene_name[$10] = $16
}

END{
for (gene in genes) {
print \
chr[gene]"\t" \
"ENSEMBL\t" \
"gene\t" \
start[gene]"\t" \
end[gene]"\t" \
".\t" \
strand[gene]"\t" \
".\t" \
"gene_id "gene" " \
"transcript_id "gene" " \
"gene_type "gene_type[gene]" " \
"gene_status \"NULL\"; " \
"gene_name "gene_name[gene]" " \
"transcript_type "gene_type[gene]" " \
"transcript_status \"NULL\"; " \
"transcript_name "gene_name[gene]" " 
}
}
' 
