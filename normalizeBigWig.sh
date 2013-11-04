#!/usr/bin/env bash

# The script:
# 1) calculates summary statistcs (coverage and mean) for a bigWig file
# 2) return another bigWig file, where the signal is normalized by the total signal*10^9

if [ $# != 3 ]; then
echo "
USAGE: $0 file.bw chrom.sizes output_dir

The script:
1) calculates summary statistics (coverage and mean) for a bigWig file
2) return another bigWig file, where the signal is normalized by the total signal*10^9
Note: the extension of the bigWig file must be .bw"
exit 1
fi

# -- Read arguments --

bigWig=$1
chromSizes=$2
outDir=$3

# -- PREPARE DIRS --

# Create a directory to store the statistics of the bigwig
totalsDir=${outDir}/totals
mkdir -p $totalsDir
# Create a directory to store the normalized files
normDir=${outDir}/norm
mkdir -p $normDir


# -- BEGIN --

# Output summary statistics for the bigWig
totalsFile=$totalsDir/$(basename $bigWig .bw).tot
/users/rg/abreschi/Documents/utils/bigWigSummaryAllChroms.sh $bigWig > $totalsFile

# Extract the total signal of the bigWig
total=$(awk 'END{print $NF}' $totalsFile)
bigWigToBedGraph $bigWig stdout | awk -v total=$total '{$4=$4/total*10^9;print}' > $(basename $bigWig .bw).bedGraph
bedGraphToBigWig $(basename $bigWig .bw).bedGraph $chromSizes $normDir/$(basename $bigWig .bw).norm.bw
rm $(basename $bigWig .bw).bedGraph

exit 0
