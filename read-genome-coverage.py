#!/usr/bin/env python

# *** Oct. 24th, 2013 ***

from argparse import ArgumentParser
from datetime import datetime
from os import path
import subprocess as sp
import multiprocessing as MP
from os import remove

# Take only the primary alignments

# Compute the read coverage for different genomic regions, given a bam file
# - exons
# - introns
# - exon-intron junctions
# - intergenic

def countReads(gtf, key):
	res = sp.Popen("samtools view -b -F 260 %s | intersectBed -abam stdin -b %s -wa -wb -f 1 -bed | wc -l" %(args.bam, gtf), shell=True, stdout=sp.PIPE)
	Q.put((key, int(res.stdout.readlines()[0].strip())))

def splitMaps(gtf, key):
	res = sp.Popen("samtools view -b -F 260 %s | bamToBed -i stdin -split | intersectBed -a stdin -b %s -f 1 -wao | groupBy -i stdin -g 4 -opCols 4,10 -ops count,min | awk '$2>1 && $3>0' | wc -l" %(args.bam, gtf), shell=True, stdout=sp.PIPE)
	Q.put((key, int(res.stdout.readlines()[0].strip())))

def countSam(key):
	res = sp.Popen("samtools view -c -F 260 %s" %(args.bam), shell=True, stdout=sp.PIPE)
	Q.put((key, int(res.stdout.readlines()[0].strip())))



# ------------------ ARGUMENT PARSING -------------------


parser = ArgumentParser(description = "Count the number of reads in genomic regions. Requires 6 CPUs")
parser.add_argument("-g", "--gtf", type=str, help="gtf with all elements (genes, transcripts and exons)")
parser.add_argument("-b", "--bam", type=str, help="bam file")
parser.add_argument("-o", "--output", type=str, help="output file name")
parser.add_argument("-I", "--ID", type=str, help="the ID of the experiment, from which the bam comes from")
#parser.add_argument("-p", "--cores", type=int, help="number of CPUs", default=1)
args = parser.parse_args()


# -------------------PRELIMINARY GTF PARSING -----------------

bn_bam = path.basename(args.bam)

print datetime.now()

## extract exons from the gtf
exon_gtf = bn_bam + "." + path.basename(args.gtf).rsplit(".", 1)[0] + ".exons.gtf"
sp.call("awk '$3 == \"exon\"' %s > %s" %(args.gtf, exon_gtf), shell=True)
print datetime.now(), "Exons extracted."

## merge the overlapping exons with mergeBed
merged_exons = path.basename(exon_gtf).rsplit(".", 1)[0] + ".merged.gtf"
sp.call("sort -k1.4,1.5n -k4n %s | mergeBed -i stdin > %s" %(exon_gtf, merged_exons), shell=True)
print datetime.now(), "Exons merged."

## extract genes from the gtf
gene_gtf = bn_bam + "."  + path.basename(args.gtf).rsplit(".", 1)[0] + ".gene.gtf"
sp.call("awk '$3 == \"gene\"' %s > %s" %(args.gtf, gene_gtf), shell=True)
print datetime.now(), "Genes extracted."

## merge the overlapping genes with mergeBed
merged_genes = path.basename(gene_gtf).rsplit(".", 1)[0] + ".merged.gtf"
sp.call("sort -k1.4,1.5n -k4n %s | mergeBed -i stdin > %s" %(gene_gtf, merged_genes), shell=True)
print datetime.now(), "Genes merged."

## extract introns with subtractBed
intron_gtf = bn_bam + "."  + path.basename(args.gtf).rsplit(".", 1)[0] + ".intron.gtf"
sp.call("subtractBed -a %s -b %s > %s" %(merged_genes, exon_gtf, intron_gtf), shell=True)
print datetime.now(), "Introns extracted."



## ------------------- SUBMIT THE READ COUNTS TO DIFFERENT PROCESSES -----------------

# Create the queue object
Q = MP.Queue()

p1 = MP.Process(target=countReads, args=(merged_exons, "exonic_reads",))
p2 = MP.Process(target=countReads, args=(intron_gtf, "intronic_reads",))
p3 = MP.Process(target=countReads, args=(merged_genes, "genic_reads",))
p4 = MP.Process(target=countSam, args=("total_reads",))
p5 = MP.Process(target=splitMaps, args=(merged_exons, "split_reads",))

for p,tag in zip((p1, p2, p3, p4, p5), ("exonic","intronic","genic","total", "split")):
	print datetime.now(), "Counting %s reads... " %tag
	p.start()

for p in (p1, p2, p3, p4, p5):
	p.join()


##------------------------ READ QUEUE RESULTS ------------------------------------------

d = {}
while not Q.empty():
	res = Q.get()
	d[res[0]] = res[1]

#d["exonic_intronic_reads"] = d["genic_reads"] - d["exonic_reads"] - d["intronic_reads"]
d["exonic_intronic_reads"] = d["genic_reads"] - d["exonic_reads"] - d["intronic_reads"] - d["split_reads"]
d["intergenic_reads"] = d["total_reads"] - d["genic_reads"]



## ------------------------ WRITE OUTPUT TO FILE ---------------------------------------


out_f = open(args.output, "w")
for k,v in d.iteritems():
	out_f.write("\t".join((args.ID, str(v), str(k))) + "\n")
out_f.close()

print datetime.now()

for f in (exon_gtf, merged_exons, gene_gtf, merged_genes, intron_gtf):
	remove(f)


exit()

