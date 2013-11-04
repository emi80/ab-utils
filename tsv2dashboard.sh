#!/bin/bash

awk 'BEGIN{FS="\t"} {if (NR==1) {for(i=1;i<=NF;i++){header[i]=$i;if($i=="FILE"){file_index=i}}} else {printf $file_index"\t"; for(i=1;i<NF;i++){if(i!=file_index){printf header[i]"="$i"; "}}printf header[i]"="$i";\n"}}' 
