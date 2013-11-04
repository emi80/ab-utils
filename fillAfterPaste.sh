#!/usr/bin/env bash

# fills up the missing values after the command paste

file=$1
fill=$2


if [[ ! "$fill" ]]; then
fill="NA"
fi

awk -v fill=$fill 'BEGIN{FS="\t"}{for(i=1;i<=NF;i++){if($i==""){$i=fill}}; print $0}' $file
