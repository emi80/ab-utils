#!/usr/bin/env python

from optparse import OptionParser
import textwrap
from collections import OrderedDict
import sys

parser = OptionParser()
parser.add_option('-f','--fasta',dest='fasta',help='fasta file. Can be stdin',metavar='file.fa')
parser.add_option('-r','--remove',dest='remove',help='list of ids',default=None,metavar='file.list')
parser.add_option('-s','--select',dest='select',help='list of ids',default=None,metavar='file.list')
options,args = parser.parse_args()

def read_fasta(fname):
	fa = open(fname) if fname != "stdin" else sys.stdin
	d = OrderedDict()
	for line in fa:
		if line.startswith('>'):
			name = line[1:].strip()
			seq = ''
			continue
		seq += line.strip()
		d[name] = seq
	return d

fasta = read_fasta(options.fasta)

if options.remove is not None:
	remove = set(line.strip() for line in open(options.remove))
	for name,seq in fasta.iteritems():
		if any(map(lambda x: x in name,remove)):
			continue
		print '>%s' %name
		print '\n'.join(textwrap.wrap(seq,70))

if options.select is not None:
	select = set(line.strip() for line in open(options.select))
	for name,seq in fasta.iteritems():
#		if any(map(lambda x: x not in name,select)):
#			continue
		if name not in select:
			continue
		print '>%s' %name
		print '\n'.join(textwrap.wrap(seq,70))
	
	
