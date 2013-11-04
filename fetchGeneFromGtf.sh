#!/usr/bin/env bash

# Given a list of gene ids, extract the gtf lines with the id in the 10th field 
# if the gene is has a dot in the gtf, that is also taken

gene=$1
gtf=$2


awk -v file=$gene 'BEGIN{while(getline<file>0) {genes[$1]}}{gene=$10;split(gene,a,"\"");split(a[2],b,".");if(a[2] in genes || b[1] in genes){print}}' $gtf
