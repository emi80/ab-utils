#!/bin/bash

awk 'BEGIN{FS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
