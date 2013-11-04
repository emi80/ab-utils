#!/bin/bash

#file=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/MakeSP/gen10.long.makeSP.exonproj.gff

if [ $# != 4 ]; then 
echo "USAGE: $0 <file> <nb_bins> <5' nb_bins> <3' nb_bins>"
exit 1
fi

file=$1
bins=$2
upstream=$3
downstream=$4

sort -k1.4n,1.5nr -k10,10 -k4n,4n -k5n,5n $file |\

awk -v file=$file -v bins=$bins -v upstream=$upstream -v downstream=$downstream '
BEGIN{
OFS="\t"
FS="\t"
while(getline<file>0){
	split($9, tags, ";")
	split(tags[1], gene_tag, " ")
	gene=gene_tag[2];
	gnlen[gene]+=$5-$4+1;
	chr[gene] = $1
	ann[gene] = $2
	elem[gene] = $3
	score1[gene] = $6
	strand[gene] = $7
	score2[gene] = $8
	info[gene] = $9
	gnstart[gene]=( (gnstart[gene]=="" || gnstart[gene]>$4) ? $4:gnstart[gene] );
	gnend[gene]=( (gnend[gene]=="" || gnend[gene]<$5) ? $5:gnend[gene] )
#	print gene, gnlen[gene], gnstart[gene], gnend[gene]
	}
prev_gene=""
}

{
split($9, tags, ";")
split(tags[1], gene_tag, " ")
gene=gene_tag[2];
exonlen=$5-$4+1

# ignore genes with cDNA shorter than the nb_bins
if (gnlen[gene] < bins) {next}

# swap up and downstream if the strand is negative
if (strand[gene] == "-") {
	tmp=downstream; downstream=upstream; upstream=tmp
	}

upext=int(gnlen[gene]*upstream/bins+0.5)         # upstream extension
downext=int(gnlen[gene]*downstream/bins+0.5)         # downstream extension

if ($4 == gnstart[gene]) {    # first exon
	# print the bins of the previous exon
	if (prev_gene != gene && prev_gene != "") {
		for (bin_id in binarray) {
			split(bin_id, bin_ids, SUBSEP)
			bin_nb = bin_ids[1]
			
			split(binarray[bin_id], bin_bp, ",")
			binstart = 1e+30; binend = -1
			for (i in bin_bp) {
				binstart = bin_bp[i]>binstart ? binstart : bin_bp[i]
				binend = bin_bp[i]<binend ? binend : bin_bp[i]
				}
	
			# revert the bin order if the strand is negative
			if (strand[prev_gene] == "-") {
				bin_nb = upstream + downstream + bins - bin_nb
				}
			
			# dont print the bin if it exceeds the number of bins due to rounding
			if (bin_nb >= upstream + downstream + bins) {
				continue
				}
			print chr[prev_gene], ann[prev_gene], elem[prev_gene], binstart, binend, score1[prev_gene], strand[prev_gene], score2[prev_gene], (info[prev_gene])" bin_number \""(bin_nb)"\";"" cDNA_length \""(gnlen[prev_gene])"\";"
			}
		}
 
	# initialize counters for the gene
	seqlen[gene] = 0 
	exon_i = 0
	offset=$4-upext+1
	split("",binarray)

	# add upstream bins (only to the first exon of course)
	for (bp=offset;bp<$4;bp++) {      # upstream extension
		bin_nb = int( (bp-offset)*bins/gnlen[gene] + 0.5 )
		binarray[bin_nb,exon_i] = ( binarray[bin_nb,exon_i]=="" ? bp : (binarray[bin_nb,exon_i])","(bp) )
		}
	}

exon_i ++

for (bp=$4;bp<=$5;bp++) {
	bin_nb = int( (bp-$4+seqlen[gene])*bins/(gnlen[gene]) + 0.5 ) + upstream 
	binarray[bin_nb,exon_i] = ( binarray[bin_nb,exon_i]=="" ? bp : (binarray[bin_nb,exon_i])","(bp) )
	}

if ($5 == gnend[gene]) {           # last exon
	
	exon_i ++

	# add downstream bins (only the last exon of course)
	for (bp=$5+1;bp<=$5+downext-1;bp++) {
		bin_nb = int( (bp-($5+1))*bins/gnlen[gene] + 0.5 ) + upstream + bins
		binarray[bin_nb, exon_i] = ( binarray[bin_nb,exon_i]=="" ? bp : (binarray[bin_nb,exon_i])","(bp) )
		}
	}	

# increment counters for next exon
seqlen[gene] += exonlen
prev_gene=gene

}

END{
for (bin_id in binarray) {
	split(bin_id, bin_ids, SUBSEP)
	bin_nb = bin_ids[1]
	
	split(binarray[bin_id], bin_bp, ",")
	binstart = 1e+30; binend = -1
	for (i in bin_bp) {
		binstart = bin_bp[i]>binstart ? binstart : bin_bp[i]
		binend = bin_bp[i]<binend ? binend : bin_bp[i]
		}

	# revert the bin order if the strand is negative
	if (strand[prev_gene] == "-") {
		bin_nb = upstream + downstream + bins - bin_nb
		}
	
	# dont print the bin if it exceeds the number of bins due to rounding
	if (bin_nb >= upstream + downstream + bins) {
		continue
		}
	print chr[prev_gene], ann[prev_gene], elem[prev_gene], binstart, binend, score1[prev_gene], strand[prev_gene], score2[prev_gene], (info[prev_gene])" bin_number \""(bin_nb)"\";"" cDNA_length \""(gnlen[prev_gene])"\";"
	}
}
'

