#!/soft/bin/python

from pylab import *
from optparse import OptionParser
from math import log, log10

parser = OptionParser()
parser.add_option('-f', '--file', dest='filename', help='', metavar='')
parser.add_option('-v', '--value', dest='value', help='the value you want to plot', metavar='')
parser.add_option('-i', '--idr', dest='idr', help='the idr threshold for different colouring', metavar='')
parser.add_option('-p', '--pseudocount', dest='psdcount', help='pseudocount for log.', metavar='', default=1)
parser.add_option('-q', '--quiet', action='store_false', dest='verbose', default=True)
options, args = parser.parse_args()


def read_gff(fname):
	l = []
	for line in open(fname):
		if line.startswith("#"):
			continue
		chr, ann, element, start, end, score1, strand, score2, mdata = line.strip().split('\t')
		mdata = dict(data.split() for data in mdata.strip(";").split("; "))
		key1, key2 = options.value+"1", options.value+"2"
		v1, v2, idr = eval(mdata[key1].strip('"')), eval(mdata[key2].strip('"')), mdata["iIDR"].strip('"')
		idr = 0 if idr == "NA" and v1 == 0 and v2 == 0 else eval(idr)
		l += [(v1, v2, idr)]
	return l

#

data = read_gff(options.filename)
p = options.psdcount

x = tuple(0.5*(log10(v1+p)+log10(v2+p)) for v1,v2,idr in data if idr > eval(options.idr))
y = tuple(log10(v1+p)-log10(v2+p) for v1,v2,idr in data if idr > eval(options.idr))

x1 = tuple(0.5*(log10(v1+p)+log10(v2+p)) for v1,v2,idr in data if idr <= eval(options.idr))
y1 = tuple(log10(v1+p)-log10(v2+p) for v1,v2,idr in data if idr <= eval(options.idr))

plot(x, y, color='black', markersize=5, linewidth=0, marker='x')
plot(x1, y1, 'o', color='red', markersize=3, linewidth=0, marker='.')


show()
