#!/usr/bin/env python

from __future__ import print_function, division, unicode_literals
import os
import sys
import json
import logging
import tempfile
import itertools
import traceback
import subprocess as sp
import multiprocessing as MP
from os.path import basename
from datetime import datetime
from argparse import ArgumentParser, FileType

PREPROC_CMDS = {
    'exon': "awk '$3 == \"exon\"' {input[0]} | sort -k1,1 -k4,4n | mergeBed -i stdin | awk 'BEGIN{{OFS=\"\t\"}}{{$(NF+1)=\"exon\";print}}' > {output}",
    'gene': "awk '$3 == \"gene\"' {input[0]} | sort -k1,1 -k4,4n | mergeBed -i stdin | awk 'BEGIN{{OFS=\"\t\"}}{{$(NF+1)=\"gene\";print}}' > {output}",
    'intron': "subtractBed -a {input[0]} -b {input[1]} | awk 'BEGIN{{OFS=\"\t\"}}{{$(NF)=\"intron\";print}}' > {output}",
    'intergenic': "complementBed -i {input[0]} -g {input[1]} | awk 'BEGIN{{OFS=\"\t\"}}{{$(NF+1)=\"intergenic\";print}}' > {output}"
}

def strfdelta(tdelta, fmt):
    d = {"days": tdelta.days}
    d["hours"], rem = divmod(tdelta.seconds, 3600)
    d["minutes"], d["seconds"] = divmod(rem, 60)
    return fmt.format(**d)

def grouper_nofill(n, iterable):
    '''list(grouper_nofill(3, 'ABCDEFG')) --> [['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
    '''
    def take():
        while 1:
            yield list(itertools.islice(iterable,n))
    return iter(take().next,[])

def preprocess(element, inputs=None):
    '''element can be one of <gene> <exon> <intron> <intergenic>'''
    log = logging.getLogger('gencov')
    #element_bed = bn_bam + "." + bn_gtf + "." + element + ".bed"
    element_bed = tempfile.mkstemp(suffix='.bed')[1]
    # if os.path.exists(element_bed):
    #     log.warning("%s bed already exists" % element.title())
    # else:
    if not inputs:
        inputs = [ args.annotation ]
    else:
        inputs = inputs[element]
    command = PREPROC_CMDS[element].format(input=inputs, output=element_bed)
    
    proc = sp.Popen(command, shell=True, stderr=sp.PIPE)
    err_msg = proc.communicate()[1]
    if err_msg:
        raise IOError(err_msg)

    log.info("%s preprocessed" % element.title())
    return element_bed

def gtf_processing(cores=1, genome=None, prefix='gencov'):
    """Annotation preprocessing. Provide a bed file with the 
    following elements:

        - projected exons
        - projected genes
        - introns
        - integenic regions

    """
    all_bed = prefix + ".all.bed"

    if not os.path.exists(all_bed) or os.stat(all_bed).st_size == 0:
        log.info("Preprocessing annotation...")
        features = ('exon', 'gene', 'intron', 'intergenic')
        pool = MP.Pool(processes=cores)
        merged_exons, merged_genes = pool.map(preprocess, features[:2])
        pool.close()
        pool.join()
        ins = {
            'intron': [merged_genes, merged_exons],
            'intergenic': [merged_genes, genome]
        }
        intron_bed, intergenic_bed = map(preprocess, features[2:], [ins, ins])

        log.info("Concatenate bed files for elements...")
        with open(all_bed, 'w') as out_bed:
            cat_all(merged_exons, merged_genes, intron_bed, intergenic_bed, out_bed=out_bed)

        for f in (merged_exons, merged_genes, intron_bed, intergenic_bed):
            os.remove(f)

    return all_bed

def cat_all(*args, **kwargs):
    out_bed = kwargs.get('out_bed', sys.stdout)
    for bed in args:    
      print(open(bed,'r').read(), end='', file=out_bed)

def get_chromosomes(genome_file):
    with open(genome_file) as genome:
        chrs = " ".join([l.split()[0] for l in genome])
    return chrs

def process_bam(bam, all_elements, chrs=None, uniq=False, all_reads=False):
    if not os.path.exists(bam):
        raise IOError("Fail to open {0!r} for reading".format(bam))
    bai = "{0}.bai".format(bam)
    if chrs and not os.path.exists(bai):
        log.info("Indexing {0}...".format(bam))
        sp.call('samtools index {0}'.format(bam), shell=True)

    log.info('Processing {0}...'.format(bam))
    command = "samtools view -u"
    if not all_reads:
        command += " -F 260"
    command += " {0}".format(bam)    
    if chrs:
        command += " {0}".format(chrs)
    if uniq:
        command = "{0} | bamtools filter -tag NH:1".format(command)
    command = "{0} | bamToBed -i stdin -cigar -bed12 | intersectBed -a stdin -b {1} -split -wao".format(command, all_elements)
    log.debug(command)
    return sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, bufsize=1)

def update_counts(element, tot_counts, cont_counts, split_counts, is_split):
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

def count_features(bed, all_elements):
    
    # Initialize
    newRead = False     # keep track of different reads
    prev_rid = None     # read id of the previous read
    is_split = False    # check if current read is a split
    element = []        # list with all elements intersecting the read
    NR = 0              # Row number (like awk)
    cont_counts = {}    # Continuous read counts
    split_counts = {}   # Split read counts
    tot_counts = {}      # Total number of reads

    # command = "intersectBed -a stdin -b %s -split -wao" % (all_elements)
    # out = sp.Popen(command, shell=True, stdout=sp.PIPE, stdin=bed.stdout, bufsize=0)

    # o = out.stdout
    o = bed.stdout

    # Iterate
    while True:
        try:
            line = o.next()
            if not line:
                continue
            NR += 1
            if 'gene' in line:
                continue
            rchr, rstart, rend, rid, rflag, rstrand, rtstart, rtend, rrgb, rbcount, rbsizes, rbstarts, achr, astart, aend, ael, covg = line.strip().split("\t")
            newRead = (rid != prev_rid)
            if (newRead) and prev_rid!=None:
                update_counts(element, tot_counts, cont_counts, split_counts, is_split)
                # Re-Initialize the counters
                element = []
            
            element.append(ael)
            prev_rid = rid
            is_split = int(rbcount) > 1
        except StopIteration:
            update_counts(element, tot_counts, cont_counts, split_counts, is_split)
            break

    return (tot_counts, cont_counts, split_counts)
    
def count_features_parallel(q,r, all_elements):
    
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

        # command = "intersectBed -a stdin -b %s -split -wao" % (all_elements)
        # out = sp.Popen(command, shell=True, stdout=sp.PIPE, stdin=sp.PIPE, bufsize=0)

        bed_lines = q.get()

        if bed_lines == "DONE":
            q.task_done()
            break

        # o = iter(out.communicate(''.join(bed_lines))[0].split('\n'))

        o = iter(bed_lines)

        # Iterate
        while True:
            try:
                line = o.next()
                if not line:
                    continue
                NR += 1
                if 'gene' in line:
                    continue
                rchr, rstart, rend, rid, rflag, rstrand, rtstart, rtend, rrgb, rbcount, rbsizes, rbstarts, achr, astart, aend, ael, covg = line.strip().split("\t")
                newRead = (rid != prev_rid)
                if (newRead) and prev_rid!=None:
                    update_counts(element, tot_counts, cont_counts, split_counts, is_split)
                    # Re-Initialize the counters
                    element = []
                
                element.append(ael)
                prev_rid = rid
                is_split = int(rbcount) > 1
            except StopIteration:
                update_counts(element, tot_counts, cont_counts, split_counts, is_split)
                break
        r.put((tot_counts, cont_counts, split_counts))
        q.task_done()
    # out.terminate()
    r.put("DONE")

def count_features_parallel2(bed):
    
    # Initialize
    newRead = False     # keep track of different reads
    prev_rid = None     # read id of the previous read
    is_split = False    # check if current read is a split
    element = []        # list with all elements intersecting the read
    NR = 0              # Row number (like awk)
    cont_counts = {}    # Continuous read counts
    split_counts = {}   # Split read counts
    tot_counts = {}      # Total number of reads

    o = iter(bed)

    # Iterate
    while True:
        try:
            line = o.next()
            if not line:
                continue
            NR += 1
            if 'gene' in line:
                continue
            rchr, rstart, rend, rid, rflag, rstrand, rtstart, rtend, rrgb, rbcount, rbsizes, rbstarts, achr, astart, aend, ael, covg = line.strip().split("\t")
            newRead = (rid != prev_rid)
            if (newRead) and prev_rid!=None:
                update_counts(element, tot_counts, cont_counts, split_counts, is_split)
                # Re-Initialize the counters
                element = []
            
            element.append(ael)
            prev_rid = rid
            is_split = int(rbcount) > 1
        except StopIteration:
            update_counts(element, tot_counts, cont_counts, split_counts, is_split)
            break
    return (tot_counts, cont_counts, split_counts)

def write_output(stats, out, output_format='tsv', json_indent=4):
    if not args.ID:
        args.ID = basename(args.bam)

    if output_format == 'tsv':
        for k, v in stats.iteritems():
            for k1, v1 in v.iteritems():
                line_array = [args.ID, k, str(k1), str(v1)]
                out.write("\t".join(line_array)+"\n")
    elif output_format == 'json':
        out.write('Total reads: %s\n' % json.dumps(stats['total'], indent=json_indent))
        out.write('Continuous reads: %s\n' % json.dumps(stats['continuous'], indent=json_indent))
        out.write('Split reads: %s\n' % json.dumps(stats['split'], indent=json_indent))

def run_parallel2(bed, chunksize=50000, cores=2):
    pool = MP.Pool(processes=cores)
    results = pool.apply(count_features_parallel2, grouper_nofill(chunksize, bed.stdout))
    pool.close()
    pool.join()
    log.info("Collect results...")
    tot = {}
    cont = {}
    split = {}
    for t,c,s in results:
        for k,v in t.items():
            tot[k] = tot.get(k,0) + v
        for k,v in c.items():
            cont[k] = cont.get(k,0) + v
        for k,v in s.items():
            split[k] = split.get(k,0) + v
    return (tot, cont, split)
    
def run_parallel(bed, all_elements, cores=2):

    q = MP.JoinableQueue()
    r = MP.Queue()

    procs = []

    for i in range(cores):
        p = MP.Process(target=count_features_parallel,args=(q, r, all_elements))
        procs.append(p)
        p.daemon = True
        p.start()
        log.debug("{0} started".format(p.name))

    log.info("Compute stats...")
    for chunk in grouper_nofill(args.chunk_size, bed.stdout):
        q.put(chunk)
    for i in range(cores):
        q.put("DONE")

    q.join()

    for p in procs:
        p.join()

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
            if proc_count == cores:
                break
            continue

        t,c,s = item
        for k,v in t.items():
            tot[k] = tot.get(k,0) + v
        for k,v in c.items():
            cont[k] = cont.get(k,0) + v
        for k,v in s.items():
            split[k] = split.get(k,0) + v

    return (tot, cont, split)

def main(args):

    bn_bam = os.path.basename(args.bam).rsplit(".", 1)[0]
    bn_gtf = os.path.basename(args.annotation).rsplit(".", 1)[0]

    start = datetime.now()
    log.info("Running with %s process%s" % (args.cores, 'es' if args.cores > 1 else ''))
    log.info("Chunk size: %s" % args.chunk_size)

    all_elements = gtf_processing(cores=args.cores, genome=args.genome, prefix=bn_bam + "." + bn_gtf)

    chrs = None if args.all_chrs else get_chromosomes(args.genome)
    bed = process_bam(args.bam, all_elements, chrs=chrs, uniq=args.uniq)

    # tot, cont, split = run_parallel2(bed, cores=args.cores, chunk_size=args.chunk_size) if args.cores > 1 else count_features(bed, all_elements)
    tot, cont, split = run_parallel2(bed, cores=args.cores, chunksize=args.chunk_size) if args.cores > 1 else count_features(bed, all_elements)

    stats_summary = {"total" : tot, "continuous" : cont, "split" : split}

    write_output(stats_summary, args.output, output_format=args.output_format)

    end = datetime.now() - start
    log.info('DONE (%s)' % strfdelta(end, "{hours}h{minutes}m{seconds}s"))

    if not args.keep:
        os.remove(all_elements)

def parse_arguments(argv):
    """ Parsing arguments """

    parser = ArgumentParser(argv, description = "Count the number of reads in genomic regions. NOTE: SAMtools and BAMtools must be installed and the BAM index must be present")
    parser.add_argument("-a", "--annotation", type=str, help="gtf with all elements (genes, transcripts and exons)", required=True)
    parser.add_argument("-g", "--genome", type=str, help="genome chromosome sizes", required=True)
    parser.add_argument("-b", "--bam", type=str, help="bam file", required=True)
    parser.add_argument("-o", "--output", type=FileType('w'), default=sys.stdout, help="output file name")
    parser.add_argument("-I", "--ID", type=str, help="the ID of the experiment, from which the bam comes from")
    parser.add_argument("-p", "--cores", type=int, help="number of CPUs", default=1)
    parser.add_argument("-r", "--records-in-ram", dest='chunk_size', type=int, help="number of records to be put in memory", default=50000)
    parser.add_argument("--keep", dest='keep', help="Do not delete the temporary files generated during the run", action='store_true', default=False)
    parser.add_argument("--uniq", dest='uniq', action='store_true', help="Only use uniquely mapped reads", default=False)
    parser.add_argument("--loglevel", dest='loglevel', help="Set the loglevel", default="info")
    parser.add_argument("--all-reads", dest='all_reads', action='store_true', help="Use all reads from the BAM file. Default: use primary alignments only ('samtools view -F 260')", default=False)
    parser.add_argument("--output-format", dest='output_format', help="Set the output format", default="tsv")
    parser.add_argument("--all-chromosomes", dest='all_chrs', action='store_true', help="Use all chromosomes from the BAM file header. Default: use only chromosomes in the genome index file.", default=False)

    return parser.parse_args()

def setup_logger():
    """ Logging setup """
    log = logging.getLogger("gencov")
    log.setLevel(logging.getLevelName(args.loglevel.upper()))
    ch = logging.StreamHandler()
    ch.setLevel = log.level
    fmt = logging.Formatter('%(asctime)s - %(message)s', '%Y-%m-%d %H:%M:%S')
    ch.setFormatter(fmt)
    log.addHandler(ch)
    return log

if __name__ == "__main__":
    """ 
    Given a bam file, compute the read coverage for different genomic regions:
        
        - exons
        - introns
        - exon-intron junctions
        - intergenic

    *** ONLY PRIMARY alignments are used ***
    """
    try:
        args = parse_arguments(sys.argv)
        log = setup_logger()
        main(args)
        exit(0)
    except Exception,err:
        log.error("Error:")
        errinfo = traceback.format_exception(sys.exc_type, sys.exc_value, sys.exc_traceback)
        log.error("".join(errinfo))
        exit(1)

