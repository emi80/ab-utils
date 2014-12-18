#!/usr/bin/env python

## LIBRARIES

import sys
from optparse import OptionParser
import gzip

## OPTION PARSING

parser = OptionParser()
parser.add_option("-I", "--input", dest="input", default="stdin",
    help="The name of the input file. The input format can be either a dashboard index file or a list of gff/gtf. leave empty for stdin.")
parser.add_option("-d", "--output_dir", dest="output_dir", help="the folder where you want to store the output file", default="./")
parser.add_option('-o','--output',dest='output',help='output file additional tags. If stdout prints on screen',default='encode',metavar='encode')
parser.add_option('-O','--organism',dest='organism',help='human or mouse',metavar='', default="human")
parser.add_option('--filters', dest='filters', help='list of variables used for filtering', default='idr')
parser.add_option('-i','--idr',dest='idr',help='select idr threshold in the interval [0,1]. A threshold of 1 means no filtering.', default="0.1", metavar='0.1')
parser.add_option('-t','--thr',dest='thr',help='lower threshold for value',default="0", metavar='0')
parser.add_option('-f','--frac_cov',dest='frac_cov',help='region covered at least',default="0", metavar='0')
parser.add_option('-F','--Frac_cov',dest='Frac_cov',help='region covered at most',default="1", metavar='1')
parser.add_option('-e','--element', dest='element',help='choose the element you want\
 to extract the matrix for: <gene>, <transcript>, <exon>, <intron>',metavar='')
parser.add_option('-s','--strip_id', dest='strip_id', action='store_true', help='strip the id removing the .[0-9]')
parser.add_option('-n','--names', default="labExpId",
    help='The metadata you want to include in the header (comma-separated). Choose "NA" for using the file name')
parser.add_option('-r','--replicates', default=False, action='store_true', help='use this if you want to use only replicates')
parser.add_option('-v','--value', default='rpkm', help='choose the values (comma-separated)\
 you want to extract: it must be a valid gff tag', metavar='')
options, args = parser.parse_args()

idr_tag = "NA" if not options.replicates else ''.join(options.idr.split('.'))
output = "%s/%s.%s.%s.%s.idr_%s.tsv" %(options.output_dir, options.output,
options.organism, options.element, options.value, idr_tag)

## FUNCTIONS

def mean(array):
    return sum(el for el in array)*1.0/len(list(array))

# ============
# Read options
# ============


##------------------------------------------------------------------------------------------------------------------------------------
## BEGIN

d = {}
cell_lines = set()
features = set()

open_input = sys.stdin if options.input == "stdin" else open(options.input)

for line in open_input:
    splitline = line.strip().strip(";").split("\t")
    if len(splitline) == 2:
        fname, meta = line.strip('; \n').split('\t')
        meta = dict(el.split('=') for el in meta.split('; '))
    elif len(splitline) == 1:
        fname = splitline[0]
        meta = {}
    else:
        print "ERROR: Wrong input format!"
        exit()
    if options.names != "NA":
        mdata_names = options.names.split(",")
        sample = '_'.join(meta.get(v, "NA") for v in mdata_names)
    else:
        sample = fname.split("/")[-1]
    cell_lines.add(sample)
    if fname.endswith('.gz'):
        my_open = gzip.open
    else:
        my_open = open
    for line in my_open(fname):
        chr, ann, segment, start, end, dot1, strand, dot2, tags = line.strip().strip(';').split('\t')
        if options.element != segment:
            features.add(segment)
            continue
        if ann == 'Cufflinks':
            continue
        if 'cuff' in tags.lower():
            continue
        tags = dict(tag.split() for tag in tags.split('; '))

        if options.element == 'transcript':
            element_id = tags['transcript_id'].strip('"')
        elif options.element == 'exon' or options.element == 'intron':
            element_id = '_'.join((chr, start, end, strand))
        else:
            element_id = tags['gene_id'].strip('"')
        if options.strip_id: element_id = element_id.split('.')[0]

        if options.replicates:
            key1=options.value+"1"
            key2=options.value+"2"
            if not tags.has_key(key1):
                cell_lines.remove(sample)
                break
            idr = tags['iIDR'].strip('"')

            if idr != 'NA' and float(idr) > float(options.idr):
                element_value = "NA"
                d.setdefault(element_id, {}).setdefault(sample,element_value)
                continue
            element_value = mean(map(lambda x: float(x.strip('"')), (tags[key1],tags[key2])))
            if element_value < float(options.thr):
                element_rpkm = "NA"
            d.setdefault(element_id, {}).setdefault(sample,element_value)
        else:
            frac_cov = float(tags.get("frac_cov","1").strip('"'))
            element_value_list = []
            for k in options.value.split(","):
                v = tags[k].strip('"')
                element_value = float(v) if v != "NA" else "NA"
                if element_value < float(options.thr) or frac_cov < float(options.frac_cov) or frac_cov > float(options.Frac_cov) :
                    element_value = 'NA'
                element_value_list += [str(element_value)]
            d.setdefault(element_id, {}). setdefault(sample, ",".join(element_value_list))


## WRITE OUTPUT

for el in features:
    sys.stderr.write("The file also contains %ss\n" %el)

f = open(output,'w') if options.output != "stdout" else sys.stdout

#f = open(output,'w')
f.write('\t'.join(sorted(cell_lines, key=lambda x: x.lower()))+'\n')
for element in sorted(d.iterkeys()):
    values = d[element]
    f.write(element+'\t')
    missing_value = ",".join(['NA']*len(options.value.split(",")))
    f.write('\t'.join(str(values.get(cell_line, missing_value)) for cell_line in sorted(cell_lines, key=lambda x: x.lower()))+'\n')
f.close()


