#!/soft/bin/python

from optparse import OptionParser
import locale


parser = OptionParser()
parser.add_option('-i', '--input', dest='input', help='the output of gt.stats', metavar='<file.stats>')
parser.add_option("-m", "--metadata", dest="metadata", help="tab separated file with metadata", metavar="<index.tsv>", default="")
parser.add_option("-f", "--fields", dest="fields", help="select the field you want in the plot title, multiple fields comma-separated", metavar="<cell>")
options, args = parser.parse_args()



#============
# Functions
#============



def read_mdata(fname):
	d = {}
	for i,line in enumerate(open(fname)):
		if i == 0:
			header = line.strip().split("\t")
			continue
		if i != 0:
			fields = line.strip().split("\t")
			k = 1 if len(fields) == len(header)+1 else 0
			id = fields[header.index("labExpId")+k]
			d[id] = "_".join(fields[header.index(field)+k] for field in options.fields.split(","))
	return d





###--------
# BEGIN
###--------


mdata = read_mdata(options.metadata) if options.metadata != "" else {}


ID = options.input.split('/')[-1].split('.')[0].split('_')[0]


for line in open(options.input):
	if line.startswith("#"):
		continue
	if line.startswith("Unmapped reads"):
		nomap = int(line.split()[-1].strip())
		continue
	if line.startswith("Multiple mapped reads"):
		mmap_reads = int(line.split()[-1].strip())
		continue
	if line.startswith("Uniquely mapped"):
		uniq_reads = int(line.split()[-1].strip())
		continue
	if line.startswith("Non-splice reads"):
		genom_reads = int(line.split()[-1].strip())
		continue
	if line.startswith("Splice reads"):
		split_reads = int(line.split()[-1].strip())
		continue
	

initial_reads = nomap + mmap_reads + uniq_reads
mapped_reads = mmap_reads + uniq_reads


l = [
("initial_reads" , initial_reads),
("mapped_reads" , mapped_reads),
("unique_reads" , uniq_reads),
("mmap_reads" , mmap_reads),
("genom_reads" , genom_reads),
("split_reads" , split_reads)
#("split1_reads" ,split1_reads),
#("split2_reads" ,split2_reads),
#("split3_reads" ,split3_reads),
#("split4p_reads" ,split4p_reads)
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
ax.set_title("%s:%s" %(ID,mdata.get(ID,"")))
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

plt.savefig(ID+"."+mdata.get(ID,"")+'.pdf', bbox_inches='tight')

#plt.show()

 




