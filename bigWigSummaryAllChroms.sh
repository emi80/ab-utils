#!/bin/bash

bigWig=$1

if [ $# != 1 ]; then
	echo "USAGE: $0 file.bigWig"
	echo "Column headers are: chr, length, covg, mean, sum"
fi

bigWigInfo $bigWig -chroms |\
awk '{if ($1~/chromCount/){out=1;N=$2;i=1;next}; if (out){print; i++;}; if(i>N){out=0;next}}' |\
while read chr index size; do
	covg=$(bigWigSummary $bigWig $chr 0 $size 1 -type=coverage)
	mean=$(bigWigSummary $bigWig $chr 0 $size 1 -type=mean)
	echo "hello" | awk -v chr=$chr -v size=$size -v covg=$covg -v mean=$mean 'BEGIN{OFS="\t"}{print chr, size, covg, mean, mean*covg*size}'
done |\
sort -k 1,1 |\
awk 'BEGIN{OFS="\t"}{print; size+=$2; covg+=$3*$2; sum+=$5}END{print "total", size, covg/size, sum/size, sum}'


exit 0
