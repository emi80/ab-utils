#!/usr/bin/env bash

f=$1
m=$2

if [[ $# != 2 ]]; then
	echo "USAGE: $0 selected_columns matrix.file.tsv"
	echo "The selected_columns are \":\" separated"
	exit 1
fi


awk -v f=$f '
BEGIN {
	FS=OFS="\t";
	split(f, a, ":");
	for (i=1;i<=length(a);i++) {
		sel[a[i]]
	}
}

NR == 1 {
	split($0, header, "\t");
	n = NF
	for (i=1;i<=length(header);i++) {
		if (header[i] in sel) {
			c++
			if (c==1) {
				printf header[i]
			} else {
				printf "\t"header[i]
			}
		}
	}
	printf "\n"
}

NR > 1 {
	c = 0
	if (NF > n+1) {
		exit 1
	}
	offset = (n == NF ? 0 : 1)
	if (offset == 1) {
		printf $1"\t"
	}
	for (i=1+offset;i<NF;i++) {
		if (header[i-offset] in sel) {
			c++
			if (c==1) {
				printf $i
			} else {
				printf "\t"$i
			}
		}
	}
	printf "\n"
}
' $m


