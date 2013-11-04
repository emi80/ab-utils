#!/bin/bash

cat - |\

awk '$3=="exon"{
print \
"chr"$1"\t" \
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
