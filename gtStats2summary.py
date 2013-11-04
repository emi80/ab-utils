#!/soft/bin/python

from optparse import OptionParser
import locale


parser = OptionParser()
parser.add_option('-f', '--file', dest='filename', help='the output of gt.stats', metavar='<file.stats>')
parser.add_option('-p', '--paired-end', action='store_true', dest='paired_end', default=False)
parser.add_option('-t', '--read-total', dest='initial_reads', help='number of initial reads', metavar='<int>')
options, args = parser.parse_args()


###
# BEGIN
###

ID = options.filename.split('/')[-1].split('.')[0].split('_')[0]

ref = ''
for line in open(options.filename):
	if "num.reads" in line.lower():
		if not options.initial_reads:
			initial_reads = int(line.strip().split()[-1])
		else:
			initial_reads = int(options.initial_reads)
		reads_in_file = int(line.strip().split()[-1])
		continue
	if "reads.mapped " in line.lower() or '--> num.mapped' in line.lower():
		mapped_reads = int(line.strip().split()[-2])
		continue
	if line.startswith('MMap.ranges'):
		ref = 'mmap'
		mmap_reads = 0
		continue
	if "[1]" in line and ref == 'mmap':
		uniq_reads = int(line.strip().split()[-2])
		continue
	if "[1]" not in line and "[0]" not in line and ref =='mmap' and line.split()[0] == '-->':
		mmap_reads += int(line.strip().split()[-2])
		continue
	if line.lower().startswith('inss.ranges'):
		ref = 'insS'
		continue
	if line.lower().startswith('uniq.ranges'):
		ref = 'uniq'
		continue
	if line.startswith('SM.Num.Splitted.Segments'):
		split_reads = int(line.strip().split()[1])
		continue
	if line.startswith("SM.Num.Junctions"):
		ref = 'smap'
		continue
	if "[1]" in line and ref == 'smap':
		split1_reads = int(line.strip().split()[-2])
		continue
	if "[2]" in line and ref == 'smap':
		split2_reads = int(line.strip().split()[-2])
		continue
	if "[3]" in line and ref == 'smap':
		split3_reads = int(line.strip().split()[-2])
		continue
	if "(3,inf)" in line and ref == 'smap':
		split4p_reads = int(line.strip().split()[-2])
		break


if options.paired_end:
	mapped_reads = mapped_reads * 2
	uniq_reads = uniq_reads * 2
	mmap_reads = mmap_reads * 2

l = [
("initial_reads" , initial_reads),
("mapped_reads" , mapped_reads),
("unique_reads" ,uniq_reads),
("mmap_reads" ,mmap_reads),
("genom_reads" , reads_in_file-split_reads),
("split_reads" ,split_reads),
("split1_reads" ,split1_reads),
("split2_reads" ,split2_reads),
("split3_reads" ,split3_reads),
("split4p_reads" ,split4p_reads)
]



#############################
# --> plotting...

import matplotlib.pyplot as plt
import numpy as np

N = len(l)

ind = np.arange(N)
width = 0.2

fig = plt.figure()
ax = fig.add_subplot(111)
rects = ax.bar(ind, tuple(n[1]*100.0/initial_reads for n in l), width, color='blue')

##add some
ax.set_ylabel('% with respect to the initial number of reads')
ax.set_title(ID)
#ax.set_xticks(ind+width)
ax.set_xticklabels( tuple(n[0] for n in l), rotation=45 )

# convert integer to string with commas
locale.setlocale(locale.LC_ALL, 'en_US.utf8')
# write labels for each bar
for rect,reads in zip(rects,l):
	read_number = locale.format('%d', reads[1], True)
	ax.text(rect.get_x()+rect.get_width()/2., 1.05*rect.get_height(), '%s'%(read_number), ha='center',va='bottom', rotation=45)


plt.ylim((0,120))
plt.xlim((-.5,N))
plt.xticks(np.arange(0, N, 1.0), ha='right')


#ax.legend( (rects1[0], rects2[0]), ('Men', 'Women') )
#
def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height), 
                ha='center', va='bottom')
#
#autolabel(rects)
#autolabel(rects2)

plt.savefig(ID+'.pdf', bbox_inches='tight')

#plt.show()

 




