#!/usr/bin/env bash

# This script extracts the 5' or 3' extreme from a gtf file and returns a gtf file

if [ $# != 2 ]; then
echo "USAGE: $0 file.gtf <5> or <3>"
exit 1
fi


gtf=$1
boundary=$2

awk -v boundary=$boundary '
BEGIN{OFS="\t";FS="\t"}
{
if (boundary == 5) {
	if ($7 == "+") {
		start = $4;
		end = $4;
		}
	if ($7 == "-") {
		start = $5;
		end = $5;
		}
	}

if (boundary == 3) {
	if ($7 == "+") {
		start = $5;
		end = $5;
		}
	if ($7 == "-") {
		start = $4;
		end = $4;
		}
	}
$4 = start
$5 = end
print
}' $gtf



#if [[ $boundary == 5 ]]; then
#awk 'BEGIN{OFS="\t";FS="\t"}{if($7=="+"){$4=$4;$5=$4}else{$4=$5;$5=$5}; print}' $gtf
#else
#awk 'BEGIN{OFS="\t";FS="\t"}{if($7=="+"){$4=$5;$5=$5}else{$4=$4;$5=$4}; print}' $gtf
#fi


