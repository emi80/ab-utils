#!/soft/bin/python

compl_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def compl_rev(seq):
	compl_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	compl_rev_strand = ''
	for nt in dna_strand:
		compl_rev_strand = compl_dict[nt] + compl_rev_strand
	return compl_rev_strand


dna_strand = raw_input('Type the dna strand here: ').upper()
dna_strand = ''.join(dna_strand.split())
print repr(dna_strand)

#rev_strand = reduce(lambda x,y: y+x, dna_strand)


compl_strand = ''
compl_rev_strand = ''
rev_strand = ''
for nt in dna_strand:
	compl_strand += compl_dict[nt]
	compl_rev_strand = compl_dict[nt] + compl_rev_strand
	rev_strand = nt + rev_strand


print 'complementary sequence ', repr(compl_strand)
print 'complentary sequence reversed ', repr(compl_rev_strand)
print 'sequence in reverse order ', repr(rev_strand)

