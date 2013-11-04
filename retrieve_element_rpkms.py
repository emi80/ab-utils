#!/usr/bin/env python

## LIBRARIES

import sys
from optparse import OptionParser
import gzip

## OPTION PARSING

parser = OptionParser()
parser.add_option("-d", "--output_dir", dest="output_dir", help="the folder where you want to store the output file", default="./")
parser.add_option('-o','--output',dest='output',help='output file additional tags',default='encode',metavar='encode')
parser.add_option('-O','--organism',dest='organism',help='human or mouse',metavar='')
parser.add_option('-i','--idr',dest='idr',help='select idr threshold in the interval [0,1]. A threshold of 1 means no filtering.', default="0.1", metavar='0.1')
parser.add_option('-t','--thr',dest='thr',help='lower threshold for value',default="0", metavar='0')
parser.add_option('-f','--frac_cov',dest='frac_cov',help='region covered at least',default="0", metavar='0')
parser.add_option('-F','--Frac_cov',dest='Frac_cov',help='region covered at most',default="1", metavar='1')
parser.add_option('-e','--element', dest='element',help='choose the element you want\
 to extract the matrix for: <gene>, <transcript>, <exon>, <intron>',metavar='')
parser.add_option('-n','--names', default=False, action='store_true', 
	help='use it if you want the metadata as column header instead of labExpId')
parser.add_option('-r','--replicates', default=False, action='store_true', help='use this if you want to use only replicates')
parser.add_option('-v','--value', default='rpkm', help='choose the values (comma-separated)\
 you want to extract: <reads>, <RPKM>, <COSI>, <incRatio>', metavar='')
options, args = parser.parse_args()

output = "%s/%s.%s.%s.%s.idr_%s.thr_%s.names_%s.tsv" %(options.output_dir, options.output,
options.organism, options.element, options.value, ''.join(options.idr.split('.')), ''.join(options.thr.split('.')), options.names)

## FUNCTIONS

def mean(array):
	return sum(el for el in array)*1.0/len(list(array))

# ============
# Read options
# ============

if options.value == 'rpkm' or options.value == 'cosi':
	options.value = options.value.upper()

##------------------------------------------------------------------------------------------------------------------------------------
## BEGIN

d = {}
cell_lines = set()
for line in sys.stdin:
	fname, meta = line.strip('; \n').split('\t')
	meta = dict(el.split('=') for el in meta.split('; '))
	if options.names:
		sample = '_'.join((options.organism,meta['cell'],meta.get('age',''),meta['lab'],meta['rnaExtract'],meta['localization'],'_'.join(meta.get('dataVersion','').split())))
	else:
		sample = meta['labExpId'].strip("\"")
	cell_lines.add(sample)
	if fname.endswith('.gz'):
		my_open = gzip.open
	else:
		my_open = open
	for line in my_open(fname):
		chr, ann, segment, start, end, dot1, strand, dot2, tags = line.strip(';\n').split('\t')
		if ann == 'Cufflinks':
			continue
		if 'cuff' in tags.lower():
			continue
		tags = dict(tag.split() for tag in tags.split('; '))
		if options.element == 'gene':
			element_id = tags['gene_id'].strip('"').split('.')[0]
		if options.element == 'transcript':
			element_id = tags['transcript_id'].strip('"').split('.')[0]
		if options.element == 'exon' or options.element == 'intron':
			element_id = '_'.join((chr, start, end, strand))

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

f = open(output,'w')
f.write('\t'.join(sorted(cell_lines))+'\n')
for element in sorted(d.iterkeys()):
	values = d[element]
	f.write(element+'\t')
	missing_value = ",".join(['NA']*len(options.value.split(",")))
	f.write('\t'.join(str(values.get(cell_line, missing_value)) for cell_line in sorted(cell_lines))+'\n')
f.close()
