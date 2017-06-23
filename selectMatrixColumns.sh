#!/usr/bin/env bash

if [[ $# != 2 ]]; then
	echo "USAGE: $0 selected_columns matrix.file.tsv"
	echo "The selected_columns are \":\" separated"
	exit 1
fi

f="$1"
m="$2"

awk -v f="$f" '
BEGIN {
	FS=OFS="\t";
	split(f, a, ":");
	for (i=1;i<=length(a);i++) {
		sel[a[i]] = i
	}
}

NR == 1 {
	n = NF
	for (i=1;i<=n;i++) {
		if ($(i) in sel) {
			indices[sel[$(i)]] = i
 		}
	}
}
{
	offset = (n == NF ? 0 : 1)
	split($0, line, FS)
	for(j in indices) {
		$(j+offset) = line[indices[j]+offset]
	} 
	NF=length(sel)+offset
}
1
' "$m"


