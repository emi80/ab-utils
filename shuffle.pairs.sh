#!/usr/bin/env bash

f=$1
n=$2

shuf $1 -n ${n} > tmp1.txt
shuf $1 -n ${n} > tmp2.txt

paste tmp1.txt tmp2.txt 

rm tmp1.txt tmp2.txt

