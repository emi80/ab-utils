#!/bin/bash



# Use this to convert a gtf to a BED6 file
if [ $# == 0 ]; then
echo "Usage: $0 file.gtf gtf_field"
echo "Default gtf_field = gene_id"
#        cat file.gtf | ./gtfgene2bed6.sh 
exit 1
fi


GTF=$1
field="gene_id"
if [[ $2 != "" ]]; then field=$2; fi

awk -v field=$field 'BEGIN{OFS=FS="\t"}
{gsub("[\";]","",$10); 
split("", info);
split($9, t, "; ");
for(i=1;i<=length(t);i++) {split(t[i], a, " "); 
#gsub(/\.[0-9]/,"",a[2]); 
gsub(/[";]/, "", a[2]); 
info[a[1]] = a[2]}
$6=($6=="."?0:$6); print $1, $4-1, $5, info[field], $6, $7}' $GTF
