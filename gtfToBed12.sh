#!/usr/bin/env bash

if [ $# != 1 ]; then
	echo "From a complete gtf file, generates a GenePred and a bed12 file."
	echo
	echo "USAGE: $0 file.gtf"
	echo
	echo "From a complete gtf file, generates a GenePred and a bed12 file."
	echo
	exit 1
fi


gtf=$1

# To create a BED12 from a GTF
# 1) convert gtf to GenePred (If an entire gtf is given it will report all the genes as errors as it
# considers only transcripts)
~abreschi/bioprogs/Jim_Kent_source_tree/gtfToGenePred $gtf -allErrors $(basename $gtf .gtf).GenePred 2> gtfToGenePred.stderr
# 2) convert GenePred to Bed
cat $(basename $gtf .gtf).GenePred | ~abreschi/bioprogs/Jim_Kent_source_tree/genePredToBed12 > $(basename $gtf .gtf).bed12

rm gtfToGenePred.stderr

