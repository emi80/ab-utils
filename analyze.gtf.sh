#!/bin/bash

gtf=$1
genome=$2
phastCons=$3

bwtool="/users/GR/mb/apohl/bin/bwtool"

if [[ $# != 3 ]]; then
	echo "
	USAGE: $0 <gtf> <genome> <phastCons>
	
	Description: create useful files from a single annotation file. 

	working dir: ./

	<gtf>: gencode format gtf file
	<genome>: single fasta file with all the chromosome sequences
	<phastCons>: single bigWig file with all the phastCons scores
	"
	exit 1;
fi

#gtf="/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/gencode.v19.annotation.gtf"
#genome="/users/rg/projects/references/Genome/Homo_sapiens.GRCh37.chromosomes.chr.M.fa"
#phastCons="/users/rg/projects/references/Genome/H.sapiens/hg19/phastCons/phastCons46way.bw/regular.chr.phastCons46way.bw"


prefix=$(basename $gtf .gtf)

# ================================ Elements =====================================================

date
printf "Making a file for each element... " > /dev/stderr

mkdir -p Elements
awk -v prefix=$prefix '$1!~/#/{print > "Elements/"prefix"."$3".gff"}' $gtf

printf "DONE\n" > /dev/stderr


# ============================== Introns ============================================================

date
printf "Adding the introns... " > /dev/stderr
 
awk '$3=="exon"' $gtf | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f /users/rg/sdjebali/Awk/make_introns.awk > Elements/${prefix}.intron.gff

printf "done\n" > /dev/stderr


# ============================= Lists ========================================================

date
printf "Writing lists of elements... " > /dev/stderr

mkdir -p Lists
awk '$3=="gene"{split($10,a,"\"");print a[2]}' $gtf | sort -u > Lists/${prefix}.gene.txt
awk '$3=="transcript"{split($12,a,"\"");print a[2]}' $gtf | sort -u > Lists/${prefix}.transcript.txt
awk '$3=="exon"{print $1"_"$4"_"$5"_"$6"_"$8}' $gtf | sort -u > Lists/${prefix}.exon.txt
awk '{print $1"_"($4-1)"_"($5+1)"_"$7}' Elements/${prefix}.intron.gff > Lists/${prefix}.txt 
awk '$3=="gene"' $gtf | ~abreschi/Documents/utils/extract.gtf.tags.sh - gene_type | grep -vw NA | sort -u > Lists/${prefix}.gene_type.txt 

printf "done\n" > /dev/stderr

# ========================= Exon projections =======================================

date
printf "Writing exon projections... " > /dev/stderr

mkdir -p MakeSP
~sdjebali/bin/make_projected_exons_per_gene.sh $gtf MakeSP/ > MakeSP/${prefix}.exonproj.gff 2> MakeSP/make_projected_exons_per_gene.err

# Divide the projected exons by gene_type
for t in `cat Lists/${prefix}.gene_type.txt`; do grep "gene_type \"$t\";" Elements/${prefix}.gene.gff | ~abreschi/Documents/utils/extract.gtf.tags.sh - gene_id | ~abreschi/Documents/utils/fetchGeneFromGtf.sh - MakeSP/${prefix}.exonproj.gff > MakeSP/${prefix}.${t}.exonproj.gff; done


# ======================== Gtf to tsv =================================================

date
printf "Writing two-column files of different gtf tags... " > /dev/stderr

# Write parsed tsv files from the gtf
mkdir -p gtf2tsv
~abreschi/Documents/utils/extract.gtf.tags.sh Elements/${prefix}.transcript.gff gene_id,transcript_id > gtf2tsv/gene.tx.tsv
~abreschi/Documents/utils/extract.gtf.tags.sh Elements/${prefix}.gene.gff gene_id,gene_name | sed '1igene_id\tgene_name' > gtf2tsv/gene.gene_name.tsv
~abreschi/Documents/utils/extract.gtf.tags.sh Elements/${prefix}.gene.gff gene_id,gene_type | sed '1igene_id\tgene_type'> gtf2tsv/gene.gene_type.tsv
awk 'BEGIN{OFS="\t"}{split($10, a, "\""); print $1"_"$4"_"$5"_"$7, a[2]}' Elements/${prefix}.exon.gff | sort -u > gtf2tsv/exon.gene_id.tsv

printf "done\n" > /dev/stderr

# ========================= Gene type =======================================

date
printf "Split file by gene_type... " > /dev/stderr

# Separate gtf by gene_type
mkdir -p gene_type
cat Lists/${prefix}.gene_type.txt | while read a; do 
grep "gene_type \"$a\";" Elements/${prefix}.gene.gff > gene_type/${prefix}.gene.${a}.gff; 
grep "gene_type \"$a\";" Elements/${prefix}.transcript.gff > gene_type/${prefix}.transcript.${a}.gff
grep "gene_type \"$a\";" Elements/${prefix}.exon.gff > gene_type/${prefix}.exon.${a}.gff
done

printf "done\n" > /dev/stderr


# ========================= Antisense ============================================

date
printf "Annotate antisense overlap... " > /dev/stderr

# Report antisense overlap
gene_type="protein_coding"
mkdir -p antisense
ls gene_type/${prefix}.gene.* | while read a; do ~/Documents/utils/classify-antisense.sh gene_type/${prefix}.gene.${gene_type}.gff $a > antisense/${gene_type}.$(basename $a | cut -f2 -d.).antisense.tsv; done
# Check which ones have exonic overlap
ls MakeSP/* | grep -v stderr | while read a; do ~/Documents/utils/classify-antisense.sh MakeSP/${prefix}.${gene_type}.exonproj.gff $a > antisense/exonic.${gene_type}.$(basename $a | cut -f1 -d.).antisense.tsv; done

printf "done\n" > /dev/stderr

# ======================= GC content ==========================================

date
printf "Computing GC content... " > /dev/stderr

mkdir -p GC-content

# -- Exons --
~abreschi/Documents/utils/gtfgene2bed6.sh MakeSP/${prefix}.exonproj.gff | fastaFromBed -bed stdin -fi $genome -fo stdout -name -s -tab | awk 'BEGIN{OFS=FS="\t"}{split(toupper($2),a,""); for (i in a) {if(a[i]=="G"||a[i]=="C"){GC[$1]++}};len[$1]+=length($2)}END{for(gene in GC){print gene, GC[gene], len[gene], GC[gene]/len[gene]*100}}' > GC-content/${prefix}.gene.exonproj.GC.tsv &

# -- All genes -- (It takes few hours...)
~abreschi/Documents/utils/gtfgene2bed6.sh Elements/${prefix}.gene.gff | fastaFromBed -bed stdin -fi $genome -fo stdout -name -s -tab | awk 'BEGIN{OFS=FS="\t"}{split(toupper($2),a,""); for (i in a) {if(a[i]=="G"||a[i]=="C"){GC[$1]++}};len[$1]+=length($2)}END{for(gene in GC){print gene, GC[gene], len[gene], GC[gene]/len[gene]*100}}' > GC-content/${prefix}.gene.GC.tsv &

# -- Introns --
#awk -v file=GC-content/${prefix}.gene.exonproj.GC.tsv 'BEGIN{OFS=FS="\t";while(getline<file>0){exlen[$1]=$3;exGC[$1]=$2}}{inGC=$2-exGC[$1]; inlen=$3-exlen[$1]; print $1, inGC, inlen, (inlen!=0?inGC/inlen:"NA")}' GC-content/${prefix}.gene.GC.tsv > GC-content/${prefix}.intronproj.GC.tsv &

printf "done\n" > /dev/stderr


# ========================== PhastCons ========================================

date
printf "Computing average exonic phastCons... " > /dev/stderr
mkdir -p phastCons

# phastCons for genes, using projected exons
~abreschi/Documents/utils/gtfgene2bed6.sh MakeSP/${prefix}.exonproj.gff | $bwtool summary -header -keep-bed -with-sum stdin $phastCons phastCons/${prefix}.exonproj.phastCons.tsv 

awk 'BEGIN{OFS=FS="\t"} $0!~/^#/ {a[$4]+=$13; b[$4]+=$7} END {for(g in a){print g, a[g]/b[g]}}' phastCons/${prefix}.exonproj.phastCons.tsv > phastCons/dmel-all-no-analysis-r5.54.exonproj.phastCons.mean.tsv 

printf "done\n" > /dev/stderr


## === GO ===
## GO terms
#mkdir -p GO
#wget ftp://ftp.flybase.net/releases/current/precomputed_files/go/gene_association.fb.gz -O GO/gene_association.fb.gz
#wget ftp://ftp.flybase.net/releases/current/precomputed_files/ontologies/gene_ontology.obo.gz -O GO/gene_ontology.obo.gz
#zcat GO/gene_ontology.obo.gz | python ~/Documents/utils/parse.obo.py -i stdin -o GO/go

## Write the annotation of the chains of genes
#mkdir -p chains
#~/Documents/utils/gene.chains.py -a ~/Documents/db/Drosophila/dmel/flybase/Element/gene.dmel-all-no-analysis-r5.54.gtf -b ~/Documents/db/Drosophila/dmel/flybase/Element/gene.dmel-all-no-analysis-r5.54.gtf > ~/Documents/db/Drosophila/dmel/flybase/chains/chains.tsv
#
## miRNA annotation
#mkdir -p miRNA
#wget ftp://mirbase.org/pub/mirbase/CURRENT/README -O miRNA/README
#wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/dme.gff3 -O miRNA/miRBase.gff3

chmod -R gu+w ./
