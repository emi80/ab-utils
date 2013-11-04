#!/usr/bin/env python

'''Calculate the profiel of a bigWig in a given region and upstream nucleotides. The internal region is binned as specified by the user. '''

import argparse 
import subprocess as SP
import itertools as IT
import multiprocessing as MP 

# == PENDING ==

# Check that all the lines are either one site or an interval
# Check that it works properly also for one site
# Check what happens if the number of bins exceed the number of nucleotides (it is dealt internally by bigWigSummary?)



# == VARIABLES ==

onesite = False
bufsize = -1
#l1, l2, l3 = [], [], []
listbuf = 1000

# == ARGS ==

parser = argparse.ArgumentParser()
parser.add_argument("-B", "--bed", type=str, help="bed file")
parser.add_argument("-w", "--bigWig", type=str, help="bigWig file")
parser.add_argument("-u", "--upstream", type=int, help="how many nucleotides upstream")
parser.add_argument("-d", "--downstream", type=int, help="how many nucleotides downstream")
parser.add_argument("-b", "--bins", type=int, help="number of bins in the range")
parser.add_argument("-p", "--processes", type=int, help="number of CPUs", default=MP.cpu_count())
args = parser.parse_args()


# == FUNCTIONS ==

def run_bigWigSummary(bigWig, chrom, start, end, dataPoints, reverse=False):
	cmnd = "bigWigSummary %s %s %s %s %s" %(bigWig, chrom, start, end, dataPoints)
	out = SP.Popen(cmnd, shell=True, bufsize=bufsize, stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE, close_fds=True).stdout.readlines()
	l = map(float, out[0].strip().replace("n/a","0").split("\t") if out else ["0"]*dataPoints)
	if reverse: l.reverse()
	return l


# == BEGIN ==

#for i,line in enumerate(open(args.bed)):
def worker(line):
	l1, l2, l3 = [], [], []
	chr, begin, end, gene, score, strand = line.strip().split("\t")
	begin, end = int(begin), int(end)
#	if end-begin == 1: 
#		onesite = True
#		if strand == "+":
#			cmnd = "bigWigSummary %s %s %s %s %s" %(args.bigWig, chr, begin-args.upstream, begin+args.downstream, args.bins)	
#			out = SP.Popen(cmnd, shell=True, bufsize=bufsize, stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE, close_fds=True).stdout.readlines()
#			l1 += [out[0].strip().replace("n/a","0").split("\t") if out else ["0"]*(args.upstream+args.downstream)]
#		if strand == "-":
#			cmnd = "bigWigSummary %s %s %s %s %s" %(args.bigWig, chr, begin-args.downstream, begin+args.upstream, args.bins)
#			out = SP.Popen(cmnd, shell=True, bufsize=bufsize, stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE, close_fds=True).stdout.readlines()
#			l1 += list(reversed(out[0].strip().replace("n/a","0").split("\t") if out else ["0"]*(args.upstream+args.downstream)))
		
#	if end-begin > 1:
#		if onesite: 
#			print "ERROR: mixed single site with interval!" 
#			exit()
#		onesite=False

	# -- UPSTREAM --
	if strand == "+":
		l1 = run_bigWigSummary(args.bigWig, chr, begin-args.upstream, begin, args.upstream)
	else:
		l1 = run_bigWigSummary(args.bigWig, chr, end, end+args.upstream, args.upstream, reverse=True)
#		if len(l1) == listbuf:
#			l1 = [[sum(map(float, i)) for i in IT.izip(*l1)]]

	# -- DOWNSTREAM --

	if strand == "+":
		l3 = run_bigWigSummary(args.bigWig, chr, end, end+args.downstream, args.downstream)
	else:
		l3 = run_bigWigSummary(args.bigWig, chr, begin-args.downstream, begin, args.downstream, reverse=True)
#		if len(l3) == listbuf:
#			l3 = [[sum(map(float, i)) for i in IT.izip(*l3)]]
	
	# -- MIDDLE --
	
#		if end-begin == 1:
#			onesite = True
#			return l1, l2 ,l3 
		
		# !!!! Check what happens if the number of bins is higher than the number of nucleotides

	if strand == "+":
		l2 = run_bigWigSummary(args.bigWig, chr, begin, end, args.bins)
	else:
		l2 = run_bigWigSummary(args.bigWig, chr, begin, end, args.bins, reverse=True)
#		if len(l2) == listbuf:
#			l2 = [[sum(map(float, i)) for i in IT.izip(*l2)]]

	return l1, l2, l3

# == END ==

#l1 = [sum(map(float, j))/(i+1) for j in IT.izip(*l1)]
#l2 = [sum(map(float, j))/(i+1) for j in IT.izip(*l2)]
#l3 = [sum(map(float, j))/(i+1) for j in IT.izip(*l3)]


# == PRINT ==

#for i,j in enumerate(l1):
#	print i-len(l1), j 

#for i,j in enumerate(l2+l3):
#	print i, j

#for i,j in enumerate(l3):
#	print "+%s" %(i+1), j
			

# == MULTITHREADING ==

pool = MP.Pool(processes=args.processes)
result = pool.imap(worker, open(args.bed, "r"))

l1, l2, l3 = [], [], []
for (en, (a, b, c)) in enumerate(result):
	l1 = [i+j for i,j in IT.izip_longest(l1,a, fillvalue=0)]
	l2 = [i+j for i,j in IT.izip_longest(l2,b, fillvalue=0)]
	l3 = [i+j for i,j in IT.izip_longest(l3,c, fillvalue=0)]

for i,j in enumerate(l1):
	print i-len(l1), j/(en+1)
for i,j in enumerate(l2+l3):
	print i, j/(en+1)


exit()





