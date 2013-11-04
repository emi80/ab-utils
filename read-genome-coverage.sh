#!/bin/bash

# *** Mar. 7th, 2013 ***

# Compute the read coverage for different genomic regions, given a list of .bam files:
# - exons
# - introns
# - exon-intron junctions
# - intergenic

gtf=$1

# ------------------ PRINT HELP -------------------

if [[ $gtf == "" ]]; then
echo "USAGE: <bam_file> | $0 <gtf_file>"
exit
fi



# -------------------PRELIMINARY GTF PARSING -----------------

date
# extract exons from the gtf
exon_gtf="$(basename $gtf .gtf).exons.gtf"
awk '$3 == "exon"' $gtf > $exon_gtf
echo $(date) Exons extracted.
# merge the overlapping exons with mergeBed
merged_exons="$(basename $exon_gtf .gtf).merged.gtf"
sort -k1.4,1.5n -k4n $exon_gtf | mergeBed -i stdin > $merged_exons
echo $(date) Exons merged.


# extract genes from the gtf
gene_gtf="$(basename $gtf .gtf).gene.gtf"
awk '$3 == "gene"' $gtf > $gene_gtf
echo $(date) Genes extracted.
# merge the overlapping genes with mergeBed
merged_genes="$(basename $gene_gtf .gtf).merged.gtf"
sort -k1.4,1.5n -k4n $gene_gtf | mergeBed -i stdin > $merged_genes
echo $(date) Genes merged.

# extract introns with subtractBed
intron_gtf="$(basename $gtf .gtf).intron.gtf"
subtractBed -a $merged_genes -b $exon_gtf > $intron_gtf
echo $(date) Introns extracted


# ------------------- WRITE THE SCRIPT TO BE SUBMITTED -----------------

# read the bam file names from stdin
while read bam_file; do

qsub_input=$(basename $bam_file).qsub
JID="$(basename $bam_file .bam)"
ID=$(echo $JID | cut -f1 -d_ | cut -f1 -d.)

rm -f ${JID}.read-genome-coverage.err
rm -f ${JID}.read-genome-coverage.out

# write the script to be submitted
echo "
#!/bin/bash
. /etc/profile

#$ -e ${JID}.read-genome-coverage.err
#$ -o ${JID}.read-genome-coverage.out

#$ -S /bin/bash
#$ -q short,rg-el6
#$ -l virtual_free=8G
#$ -m e
#$ -M ale.breschi@gmail.com
#$ -cwd

echo 'Running on' \$HOSTNAME
date 
echo 'Creating bam with only primary alignments'
#260=256(not primary alignment)+4(unmapped)
primary_bam=\$TMPDIR/$(basename $bam_file .bam).primary.bam
samtools view -b -F 260 $bam_file > \$primary_bam
date
echo 'Running intersectBed on exons'
exonic_reads=\$(intersectBed -abam \$primary_bam -b $merged_exons -wa -wb -f 1 -bed | wc -l &)
echo 'Running intersectBed on introns'
intronic_reads=\$(intersectBed -abam \$primary_bam -b $intron_gtf -wa -wb -f 1 -bed | wc -l &)
echo 'Running intersectBed on genes'
genic_reads=\$(intersectBed -abam \$primary_bam -b $merged_genes -wa -wb -f 1 -bed | wc -l &)
echo 'Calculating total reads' 
total_reads=\$(samtools view -c \$primary_bam &)

until [ \$exonic_reads != \"\" ] && [ \$intronic_reads != \"\" ] && [ \$genic_reads != \"\" ] && [ \$total_reads != \"\" ]; do 
sleep 60s
done

exonic_intronic_reads=\$(echo  \" \$genic_reads - \$exonic_reads - \$intronic_reads \" | bc)
intergenic_reads=\$(echo \" \$total_reads - \$genic_reads \" | bc)

echo -ne \"$ID\t\${exonic_reads}\texonic\n
$ID\t\${intronic_reads}\tintronic\n
$ID\t\${exonic_intronic_reads}\texonic-intronic\n
$ID\t\${genic_reads}\tgenic\n
$ID\t\${intergenic_reads}\tintergenic\n
$ID\t\${total_reads}\ttotal\n
\" > ${JID}.read-genome-coverage.stats

date

" > $qsub_input

# launch it on the cluster
qsub -N $ID $qsub_input

done



