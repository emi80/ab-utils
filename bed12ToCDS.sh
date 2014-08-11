#!/usr/bin/env bash

bed12=$1

if [[ $# <1 ]]; then 
	echo "USAGE: $0 <bed12>"
	exit 1
fi;


# Taken from
# http://onetipperday.blogspot.com.es/2012/11/get-intron-utr-cds-from-bed12-format.html


awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; if($7==$8) next; j=0; for(i=1;i<length(a);i++) if(($2+b[i]+a[i])>$7 && ($2+b[i])<$8) {j++; start=$2+b[i]-$7; size=a[i]; if(($2+b[i])<=$7) {start=0;size=size-($7-($2+b[i]));} if(($2+a[i]+b[i])>=$8) {size=size-($2+a[i]+b[i]-$8);} A=A""size",";B=B""start",";} print $1,$7,$8,$4,$5,$6,$7,$8,$9,j,A,B;}' $bed12
