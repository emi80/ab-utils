#!/usr/bin/env bash

wigFix=$1


awk 'BEGIN{OFS="\t"}
{
if($0~/fixedStep/){
split($2,a,"=");split($3,b,"=");split($4,c,"=");chr=a[2];start=b[2];step=c[2]
}else{
print chr,start-1,start,$1;start+=step
}
}' $file
