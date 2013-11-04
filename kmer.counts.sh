#!/usr/bin/env bash

fasta=$1
k=$2


if [ $k == 1 ]; then
echo -e {A,C,G,T} | tr ' ' "\n" > kmers.txt
fi

if [ $k == 2 ]; then
echo -e {A,C,G,T}{A,C,G,T} | tr ' ' "\n" > kmers.txt
fi

if [ $k == 3 ]; then
echo -e {A,C,G,T}{A,C,G,T}{A,C,G,T} | tr ' ' "\n" > kmers.txt
fi

if [ $k == 4 ]; then
echo -e {A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T} | tr ' ' "\n" > kmers.txt
fi

if [ $k == 5 ]; then
echo -e {A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T} | tr ' ' "\n" > kmers.txt
fi

if [ $k == 6 ]; then
echo -e {A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T} | tr ' ' "\n" > kmers.txt
fi



awk -v file=kmers.txt -v k=$k 'BEGIN{OFS="\t";while(getline<file>0){kmers[$1]}}
{if ($0~/>/) {gsub(/>/,"",$1); name=$1}
else{c+=length($1);for(kmer in kmers){pat[kmer, name]=0}; for(i=1;i<=length($1)-7;i++){pat[substr($1,i,k), name]++}}}
END{
for(p in pat){
split(p,a,SUBSEP)
print a[1],a[2],pat[p]}}' $fasta

rm kmers.txt


