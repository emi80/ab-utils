#!/usr/bin/env bash

if [[ $# != 2 ]]; then
echo "USAGE: $0 <bwtoolMatrix> <offset>"
exit 1
fi


bwtoolMatrix=$1
offset=$2



awk -v offset=$offset 'BEGIN{OFS="\t"}{for(i=1;i<=NF;i++){c[i]+=($i=="-nan"?0:$i)}}END{for(i=1;i<=length(c);i++){print i-offset,c[i]/NR}}' $bwtoolMatrix
