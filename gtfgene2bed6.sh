#!/bin/bash

# Use this to convert a gtf to a BED6 file
if [ $# == 0 ]; then
echo  "Usage: $0 file.gtf"
#        cat file.gtf | ./gtfgene2bed6.sh 
exit 1
fi


GTF=$1

awk '{gsub("[\";]","",$10); gsub(/\.[0-9]+/,"",$10); $6=($6=="."?0:$6); printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1,$4-1,$5,$10,$6,$7}' $GTF
