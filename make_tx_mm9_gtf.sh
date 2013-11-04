#!/bin/bash

cat - |\

awk '
{
genes[$12] = $12
chr[$12] = "chr"$1
if (start[$12] == "" || start[$12] > $4 ) {start[$12] = $4}
if (end[$12] == "" || end[$12] < $5) {end[$12] = $5}
strand[$12] = $7
gene_id[$12] = $10
gene_type[$12] = $18
gene_name[$12] = $16
transcript_type[$12] = $2
transcript_name[$12] = $20
}

END{
for (gene in genes) {
print \
chr[gene]"\t" \
"ENSEMBL\t" \
"transcript\t" \
start[gene]"\t" \
end[gene]"\t" \
".\t" \
strand[gene]"\t" \
".\t" \
"gene_id "gene_id[gene]" " \
"transcript_id "gene" " \
"gene_type "gene_type[gene]" " \
"gene_status \"NULL\"; " \
"gene_name "gene_name[gene]" " \
"transcript_type \""transcript_type[gene]"\"; " \
"transcript_status \"NULL\"; " \
"transcript_name "transcript_name[gene]" " 
}
}
' 
