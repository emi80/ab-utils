#!/bin/bash

gtf=$1
wigFix=$2

if [ $# != 2 ]; then
echo "USAGE: $0 <file.gtf> <file.wigFix>"
exit 1;
fi

awk 'BEGIN{OFS="\t"}
{
if($0~/fixedStep/){
split($2,a,"=");split($3,b,"=");split($4,c,"=");chr=a[2];start=b[2];step=c[2]
}else{
print chr,start-1,start,$1;start+=step
}
}' $wigFix |\
intersectBed -a stdin -b $gtf -wa -wb |\
awk '{gene=$14; gsub("[\";]","",gene); gsub("\\.[0-9]+","",gene)
score[gene]+=$4;n[gene]++}
END{for(gene in score){print gene,score[gene],n[gene],score[gene]/n[gene]}}'
