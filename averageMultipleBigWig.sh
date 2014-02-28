#!/bin/bash

out=$1
count=0
chromsizes=$2

if [ $# != 2 ]; then
	echo "USAGE: cat file.list.txt | $0 out_prefix chrom.sizes"
	exit 1
fi

while read f; do
	(( count++ ))
	fileStr="$f $fileStr"
done


/usr/bin/time -o ${out}.time bigWigMerge $fileStr stdout 2>/dev/null | awk -v T=$count '{$4=$4/T; print}' > ${out}.bedGraph
/software/rg/el6.3/bin/bedGraphToBigWig ${out}.bedGraph $chromsizes ${out}.bw
rm ${out}.bedGraph

exit 0


