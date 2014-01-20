#!/usr/bin/env python

# *** Oct. 24th, 2013 ***
from argparse import ArgumentParser, FileType
from datetime import datetime
import subprocess as sp
import multiprocessing as MP
import sys
from os import path, remove

# Take only the primary alignments

# Compute the read coverage for different genomic regions, given a bam file
# - exons
# - introns
# - exon-intron junctions
# - intergenic
all_bed = ""

def grouper_nofill(n, iterable):
    '''list(grouper_nofill(3, 'ABCDEFG')) --> [['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
    '''
    it=iter(iterable)
    def take():
        while 1:
            yield [r for r in list(itertools.islice(it,n))]
    return iter(take().next,[])


def project(element):
    '''element can be one of <gene> <exon>'''
    merged_element = bn_bam + "." + bn_gtf + "." + element + ".merged.bed"
    if path.exists(merged_element):
        print element, "merged bed already exists"
    else:
        sp.call("awk '$3 == \"%s\"' %s | sort -k1,1 -k4,4n | mergeBed -i stdin | awk 'BEGIN{OFS=\"\t\"}{$(NF+1)=\"%s\";print}' > %s" %(element, args.annotation, element, merged_element), shell=True)
        print datetime.now(), element, "merged."
    return merged_element

def count_features(q, r):
    while True:
        # Initialize
        newRead = False     # keep track of different reads
        prev_rid = None     # read id of the previous read
        is_split = False    # check if current read is a split
        element = []        # list with all elements intersecting the read
        NR = 0              # Row number (like awk)
        cont_counts = {}    # Continuous read counts
        split_counts = {}   # Split read counts
        tot_counts = {}      # Total number of reads


        command = "intersectBed -a stdin -b %s -split -wao" % (all_bed)
        out = sp.Popen(command, shell=True, stdout=sp.PIPE, stdin=sp.PIPE)

        bed_lines = q.get()
        o = iter(out.communicate(''.join(bed_lines))[0].split('\n'))

        # Iterate
        while True:
            line = o.next()
            NR += 1
            if 'gene' in line:
                continue
            if line != "":
                rchr, rstart, rend, rid, rflag, rstrand, rtstart, rtend, rrgb, rbcount, rbsizes, rbstarts, achr, astart, aend, ael,covg = line.strip().split("\t")
            newRead = (rid != prev_rid)
            if (newRead or line=="") and prev_rid!=None:
                elem='total'
                tot_counts[elem] = tot_counts.get(elem,0) + 1
                if is_split:
                    split_counts['total'] = split_counts.get('total',0) + 1
                    if len(element) > 1:
                        if len(set(element)) == 1:
                            elem = element[0]
                        else:
                            if 'intergenic' in element:
                                elem = 'others'
                            else:
                                elem = 'exonic_intronic'
                    else:
                        elem = element[0]

                    split_counts[elem] = split_counts.get(elem, 0) + 1

                else:
                    cont_counts['total'] = cont_counts.get('total', 0) + 1
                    if len(element) > 1:
                        if 'intergenic' in element:
                            elem = 'others'
                        else:
                            elem = 'exonic_intronic'
                    else:
                        elem = element[0]

                    cont_counts[elem] = cont_counts.get(elem, 0) + 1

                # Re-Initialize the counters
                element = []

                if line == "":
                    break
            if line != "":
                element.append(ael)
                prev_rid = rid
                is_split = int(rbcount) > 1
        #   break
        r.put((tot_counts, cont_counts, split_counts))
        q.task_done()

if __name__ == "__main__":
    # ------------------ ARGUMENT PARSING -------------------


    parser = ArgumentParser(description = "Count the number of reads in genomic regions. Requires 6 CPUs")
    parser.add_argument("-a", "--annotation", type=str, help="gtf with all elements (genes, transcripts and exons)")
    parser.add_argument("-g", "--genome", type=str, help="genome fasta")
    parser.add_argument("-b", "--bam", type=str, help="bam file")
    parser.add_argument("-o", "--output", type=FileType('w'), default=sys.stdout, help="output file name")
    parser.add_argument("-I", "--ID", type=str, help="the ID of the experiment, from which the bam comes from")
    parser.add_argument("-p", "--cores", type=int, help="number of CPUs", default=1)
    parser.add_argument("-r", "--records-in-ram", dest='chunk_size', type=int, help="number of records to be put in memory", default=10000)
    args = parser.parse_args()


    # -------------------PRELIMINARY GTF PARSING -----------------

    bn_bam = path.basename(args.bam).rsplit(".", 1)[0]
    bn_gtf = path.basename(args.annotation).rsplit(".", 1)[0]


    print datetime.now()

    # Exon projections
    merged_exons = project("exon")

    # Gene projections
    merged_genes = project("gene")

    # Introns
    intron_bed = bn_bam + "." + bn_gtf + ".intron.bed"
    if path.exists(intron_bed):
        print "Intron bed already exists"
    else:
        sp.call("subtractBed -a %s -b %s | awk 'BEGIN{OFS=\"\t\"}{$(NF)=\"intron\";print}' > %s" %(merged_genes, merged_exons, intron_bed), shell=True)
        print datetime.now(), "intron extracted."

    # Intergenic
    intergenic_bed = bn_bam + "." + bn_gtf + ".intergenic.bed"
    if path.exists(intergenic_bed):
        print "Intergenic bed already exists"
    else:
        sp.call("complementBed -i %s -g %s | awk 'BEGIN{OFS=\"\t\"}{$(NF+1)=\"intergenic\";print}' > %s" %(merged_genes, args.genome, intergenic_bed), shell=True)
        print datetime.now(), "intergenic extracted"

    all_bed = bn_bam + "." + bn_gtf + ".all.bed"
    sp.call("cat %s %s %s %s > %s" %(merged_exons, merged_genes, intron_bed, intergenic_bed, all_bed), shell=True)
    print datetime.now(), "cat all bed"


    ## -----------------------------------------------------------------------------------

    bed = sp.Popen("samtools view -b -F 260 %s | bamToBed -i stdin -cigar -bed12" %(args.bam), shell=True, stdout=sp.PIPE)

    import itertools
    import sys
    print "Running with %s processes" % args.cores
    print "Chunk size:", args.chunk_size

    q = MP.JoinableQueue()
    r = MP.Queue()

    for i in range(args.cores):
        p = MP.Process(target=count_features,args=(q, r))
        p.daemon = True
        #print "Starting", p.name
        p.start()

    for chunk in grouper_nofill(args.chunk_size, bed.stdout):
        q.put(chunk)

    q.join()

    ## ------------------------ GATHER RESULTS ---------------------------------------

    tot = {}
    cont = {}
    split = {}

    while not r.empty():
        t, c, s = r.get()
        for k,v in t.items():
            tot[k] = tot.get(k,0) + v
        for k,v in c.items():
            cont[k] = cont.get(k,0) + v
        for k,v in s.items():
            split[k] = split.get(k,0) + v


    ## ------------------------ WRITE OUTPUT TO FILE ---------------------------------------

    import json
    out = args.output
    out.write('Total reads: %s\n' % json.dumps(tot, indent=4))
    out.write('Continuous reads: %s\n' % json.dumps(cont, indent=4))
    out.write('Split reads: %s\n' % json.dumps(split, indent=4))

    print datetime.now()
    print 'DONE'

    for f in (merged_exons, merged_genes, intron_bed, intergenic_bed, all_bed):
        remove(f)

    exit()

