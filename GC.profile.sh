#!/usr/bin/env bash

echo "
Output a matrix of 1 and 0 at each nucleotide, if it is G/C or A/T, respectively. 
Take as input a single nucleotide bed and extend of the specified interval
" > /dev/stderr


bed=$1
fasta=$2
left=$3
right=$4

if [[ $# < 2 ]]; then
	echo "USAGE: $0 <file.bed> <file.fa> left right"
	exit 1
fi


### BEGIN ###

cat $bed |\
awk -v left=$left -v right=$right 'BEGIN{OFS=FS="\t"}
{start = ($6 == "+" ? $2-left : $2-right);
end = ($6 == "+" ? $3+right : $3+left);
print $1, start, end, $4, $5, $6}' |\
fastaFromBed -bed $bed -fi $fasta -fo stdout -name -s |\
awk '$0!~/^>/{
split($1, a, "");
for(i=1;i<length(a);i++){
out = (a[i]=="C"||a[i]=="G" ? 1 : 0);
printf out"\t"
}
out = (a[i]=="C"||a[i]=="G" ? 1 : 0);
printf out"\n"
}' 


