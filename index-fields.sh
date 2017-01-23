#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
