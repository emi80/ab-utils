#!/bin/bash

cat - |\

awk '$3=="exon"{
if ($1 !~/^Contig/ && $1 !~/^Ultra/) {
	chr[$12] = "chr"$1
} else {chr[$12] = $1}
print \
chr[$12]"\t" \
"ENSEMBL\t" \
"exon\t" \
$4"\t" \
$5"\t" \
".\t" \
$7"\t" \
".\t" \
"gene_id "$10" " \
"transcript_id "$12" " \
"gene_type "$18" " \
"gene_status \"NULL\"; " \
"gene_name "$16" " \
"transcript_type \""$2"\"; " \
"transcript_status \"NULL\"; " \
"transcript_name "$20" " 
}

' 
