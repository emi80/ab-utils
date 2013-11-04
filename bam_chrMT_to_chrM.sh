#!/bin/bash

# *** Mar. 9th, 2013 ***

# replace chrMT with chrM at every occurrency in a bam



# ------------------- WRITE THE SCRIPT TO BE SUBMITTED -----------------

# read the bam file names from stdin
while read bam_file; do

qsub_input=$(basename $bam_file).qsub
JID="$(basename $bam_file .bam)_chrMT"
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
#$ -q rg-el6,default
#$ -l virtual_free=8G
#$ -m e
#$ -M ale.breschi@gmail.com
#$ -cwd

echo 'Running on' \$HOSTNAME
date 

samtools view -h $bam_file | sed 's/chrMT/chrM/g' | samtools view -bS - > $(basename $bam_file .bam).chrM.bam

date

" > $qsub_input

# launch it on the cluster
qsub -N $ID $qsub_input

done



