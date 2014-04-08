#!/usr/bin/env python

# *** Feb. 28th, 2014 ***
from argparse import ArgumentParser, FileType
from datetime import datetime
import subprocess as sp
import multiprocessing as MP
import sys
import os
import traceback

# Take only the primary alignments

# Compute the read coverage for different genomic regions, given a bam file
# - exons
# - introns
# - exon-intron junctions
# - intergenic
all_bed = ""

def strfdelta(tdelta, fmt):
    d = {"days": tdelta.days}
    d["hours"], rem = divmod(tdelta.seconds, 3600)
    d["minutes"], d["seconds"] = divmod(rem, 60)
    return fmt.format(**d)

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
    import logging
    log = logging.getLogger('gencov')
    merged_element = bn_bam + "." + bn_gtf + "." + element + ".merged.bed"
    if os.path.exists(merged_element):
        log.warning("%s merged bed already exists" % element.title())
    else:
        sp.call("awk '$3 == \"%s\"' %s | sort -k1,1 -k4,4n | mergeBed -i stdin | awk 'BEGIN{OFS=\"\t\"}{$(NF+1)=\"%s\";print}' > %s" %(element, args.annotation, element, merged_element), shell=True)
        log.info("%s merged" % element.title())
    return merged_element

def count_features(q,r):
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
        out = sp.Popen(command, shell=True, stdout=sp.PIPE, stdin=sp.PIPE, bufsize=0)

        bed_lines = q.get()

        if bed_lines == "DONE":
            q.task_done()
            break

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
    r.put("DONE")

if __name__ == "__main__":
    # ------------------ ARGUMENT PARSING -------------------

    import itertools
    import sys
    import logging

    parser = ArgumentParser(description = "Count the number of reads in genomic regions. NOTE: SAMtools and BAMtools must be installed and the BAM index must be present")
    parser.add_argument("-a", "--annotation", type=str, help="gtf with all elements (genes, transcripts and exons)")
    parser.add_argument("-g", "--genome", type=str, help="genome chromosome sizes")
    parser.add_argument("-b", "--bam", type=str, help="bam file")
    parser.add_argument("-o", "--output", type=FileType('w'), default=sys.stdout, help="output file name")
    parser.add_argument("-I", "--ID", type=str, help="the ID of the experiment, from which the bam comes from")
    parser.add_argument("-p", "--cores", type=int, help="number of CPUs", default=1)
    parser.add_argument("-r", "--records-in-ram", dest='chunk_size', type=int, help="number of records to be put in memory", default=50000)
    parser.add_argument("--loglevel", dest='loglevel', help="Set the loglevel", default="info")
    parser.add_argument("--uniq", dest='uniq', action='store_true', help="Only use uniquely mapped reads", default=False)
    parser.add_argument("--keep", dest='keep', help="Do not delete the temporary files generated during the run", action='store_true', default=False)

    args = parser.parse_args()

    # set up logger
    log = logging.getLogger("gencov")
    log.setLevel(logging.getLevelName(args.loglevel.upper()))
    ch = logging.StreamHandler()
    ch.setLevel = log.level
    fmt = logging.Formatter('%(asctime)s - %(message)s', '%Y-%m-%d %H:%M:%S')
    ch.setFormatter(fmt)
    log.addHandler(ch)

    bn_bam = os.path.basename(args.bam).rsplit(".", 1)[0]
    bn_gtf = os.path.basename(args.annotation).rsplit(".", 1)[0]

    start = datetime.now()
    log.info("Running with %s process%s" % (args.cores, 'es' if args.cores > 1 else ''))
    log.info("Chunk size: %s" % args.chunk_size)

    # -------------------PRELIMINARY GTF PARSING -----------------
    log.info("Parsing annotation...")

    pool = MP.Pool(processes=args.cores)
    # Exon and gene projections
    features = ('exon', 'gene')
    merged_exons, merged_genes = pool.map(project, features)
    pool.close()
    pool.join()

    # Introns
    intron_bed = bn_bam + "." + bn_gtf + ".intron.bed"
    if os.path.exists(intron_bed):
        log.warn( "Intron bed already exists")
    else:
        sp.call("subtractBed -a %s -b %s | awk 'BEGIN{OFS=\"\t\"}{$(NF)=\"intron\";print}' > %s" %(merged_genes, merged_exons, intron_bed), shell=True)
        log.info("Intron extracted")

    # Intergenic
    intergenic_bed = bn_bam + "." + bn_gtf + ".intergenic.bed"
    if os.path.exists(intergenic_bed):
        log.warn("Intergenic bed already exists")
    else:
        sp.call("complementBed -i %s -g %s | awk 'BEGIN{OFS=\"\t\"}{$(NF+1)=\"intergenic\";print}' > %s" %(merged_genes, args.genome, intergenic_bed), shell=True)
        log.info("Intergenic extracted")

    all_bed = bn_bam + "." + bn_gtf + ".all.bed"
    sp.call("cat %s %s %s %s > %s" %(merged_exons, merged_genes, intron_bed, intergenic_bed, all_bed), shell=True)
    log.info("Cat all bed files...")

    ## -----------------------------------------------------------------------------------

    bam = args.bam
    chrs = " ".join(sorted([ l.split()[0] for l in open(args.genome) ]))

    if not os.path.exists("{0}.bai".format(args.bam)):
        try:
            tmpbam = os.path.join(os.environ.get("TMPDIR","/tmp"), os.path.basename(bam))
            if not os.path.exists(tmpbam):
                log.debug("Create symlink to {0} in {1}".format(bam, os.path.dirname(tmpbam)))
                os.symlink(os.path.abspath(bam), tmpbam)
            bam = tmpbam
        except Exception,err:
            log.error("Error indexing the bamfile:")
            errinfo = traceback.format_exception(sys.exc_type, sys.exc_value, sys.exc_traceback)
            log.error("".join(errinfo))
        else:
            log.debug("Index {0}".format(bam))
            sp.call('samtools index {0}'.format(bam), shell=True)

    log.debug('Reading {0}'.format(bam))
    command = "samtools view -b -F 260 {0} {1}".format(bam, chrs)
    if args.uniq:
        command = "{0} | bamtools filter -tag NH:1".format(command)
    command = "{0} | bamToBed -i stdin -cigar -bed12".format(command)
    log.debug(command)
    bed = sp.Popen(command, shell=True, stdout=sp.PIPE, bufsize=1)

    q = MP.JoinableQueue()
    r = MP.Queue()

    procs = []

    for i in range(args.cores):
        p = MP.Process(target=count_features,args=(q, r))
        procs.append(p)
        p.daemon = True
        p.start()
        log.debug("{0} started".format(p.name))

    log.info("Compute stats...")
    for chunk in grouper_nofill(args.chunk_size, bed.stdout):
        q.put(chunk)
    for i in range(args.cores):
        q.put("DONE")

    q.join()


    ## ------------------------ COLLECT RESULTS ---------------------------------------
    log.info("Collect results...")

    tot = {}
    cont = {}
    split = {}

    proc_count = 0

    while True:
        item = r.get()

        # check if all processes finished
        if item == "DONE":
            proc_count += 1
            if proc_count == args.cores:
                break
            continue

        t,c,s = item
        for k,v in t.items():
            tot[k] = tot.get(k,0) + v
        for k,v in c.items():
            cont[k] = cont.get(k,0) + v
        for k,v in s.items():
            split[k] = split.get(k,0) + v


    ## ------------------------ WRITE OUTPUT TO FILE ---------------------------------------

    out = args.output

    summary_d = {"total" : tot, "continuous" : cont, "split" : split}
    if not args.ID:
        from os.path import basename
        args.ID = basename(args.bam)
    for k, v in summary_d.iteritems():
        for k1, v1 in v.iteritems():
            line_array = [args.ID, k, str(k1), str(v1)]
            out.write("\t".join(line_array)+"\n")

    out.close()
#    import json
#    out.write('Total reads: %s\n' % json.dumps(tot, indent=4))
#    out.write('Continuous reads: %s\n' % json.dumps(cont, indent=4))
#    out.write('Split reads: %s\n' % json.dumps(split, indent=4))

    end = datetime.now() - start
    log.info('DONE (%s)' % strfdelta(end, "{hours}h{minutes}m{seconds}s"))

    if not args.keep:
        for f in (merged_exons, merged_genes, intron_bed, intergenic_bed, all_bed):
            os.remove(f)

    exit()

