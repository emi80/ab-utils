#!/bin/bash

file=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/MakeSP/gen10.long.makeSP.exonproj.gff

sort -nk1.4,1.5r -k9,9 -k4,4n -k5,5n $file |\

awk -v file=$file -v bins=100 -v upstream=10 -v downstream=10 '
BEGIN{
OFS="\t"
FS="\t"
while(getline<file>0){
	split($9, tags, ";")
	split(tags[1], gene_tag, " ")
	gene=gene_tag[2];
	gnlen[gene]+=$5-$4+1;
	gnstart[gene]=( (gnstart[gene]=="" || gnstart[gene]>$4) ? $4:gnstart[gene] );
	gnend[gene]=( (gnend[gene]=="" || gnend[gene]<$5) ? $5:gnend[gene] )
#	print gene, gnlen[gene], gnstart[gene], gnend[gene]
	}
}

{
split($9, tags, ";")
split(tags[1], gene_tag, " ")
gene=gene_tag[2];
binlen=int(gnlen[gene]/bins+0.5)
exonlen=$5-$4+1
exonbins=int(exonlen/binlen)

if ($4 == gnstart[gene]) {    # first exon
	exonbins=int(exonlen/binlen)
	mod=exonlen%binlen      # Module of the bin
	for (i=1;i<=upstream;i++) {
		binstart=gnstart[gene]-binlen*(upstream-i+1)
		binend=binstart+binlen-1
		print $1,$2,$3,binstart,binend,$6,$7,$8,($9)" bin_number \""(i)"\";"
		}
	for (i=upstream+1;i<=exonbins+upstream+1;i++){
		binstart=gnstart[gene]+binlen*(i-upstream-1)
		binend=binstart+binlen-1
		if (i==exonbins+upstream+1 && mod!=0) {      # end of the exon with leftover bases
			binend=$5
			}
		print $1,$2,$3,binstart,binend,$6,$7,$8,($9)" bin_number \""(i)"\";"
		}
	print "MOD", mod, i
	}

if ($4 != gnstart[gene]) {     # not initial exon
	offset=i-1
	diff=binlen-mod
	exonbins=int((exonlen-diff)/binlen)
	if (mod!=0) {
		binstart=$4
		binend=binstart+(diff-1)
		print $1,$2,$3,binstart,binend,$6,$7,$8,($9)" bin_number \""(offset)"\";"
		}
	mod=(exonlen-(diff))%binlen     # new module of the bin (DONT MOVE THIS LINE)
	for (i=offset+1;i<=exonbins+offset+1;i++) {
		binstart=(diff)+$4+binlen*(i-offset-1)
		binend=binstart+binlen-1
		if (i==exonbins+offset+1 && mod!=0) {      # end of the exon with leftover bases
		    binend=$5
			}
		print $1,$2,$3,binstart,binend,$6,$7,$8,($9)" bin_number \""(i)"\";"
		}
	print "MOD", mod, i
	diff=binlen-mod                      # new diff (DONT MOVE THIS LINE)
	}



#if ($5 == gnend[gene]) {    # terminal exon
#	if (diff!=0) {
#		
#		}
#	}
}
'



