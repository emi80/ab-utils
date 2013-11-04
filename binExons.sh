#!/bin/bash

bins=$2
file=$1

if [ $# != 2 ]; then
echo "USAGE: $0 <exons.bed> <n_bins>"
exit 1
fi



# Bin with respect to gene length
# Bin to respect to cDNA length?


awk -v bins=$bins '{
chr=$1; ex_start=$2; ex_end=$3; gene_start=$4; gene_end=$5; strand=$6
bin1 = (ex_start - gene_start)/(gene_end - gene_start)*bins 
bin2 = (ex_end - gene_start)/(gene_end - gene_start)*bins
if (int(bin1) != int(bin2)) {
	bin = bin1-int(bin1) < bin2-int(bin2) ? int(bin1) : int(bin2)
	} else { bin = int(bin1) }
if (strand == "-") {
	bin = bins -1 - bin
	}
#print bin, $0
print chr"_"ex_start+1"_"ex_end"_"strand"\t"bin
}' $file
