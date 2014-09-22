#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file")
parser.add_argument("-o", "--output", type=str, help="output file name")
args = parser.parse_args()



transcript_types = [
	"mRNA",
	"ncRNA",
	"pseudogene",
	"rRNA",
	"pre_miRNA",
	"snRNA",
	"snoRNA",
	"tRNA"
	]

valid_chr = [
	"chr2L",
	"chr2R",
	"chr3L",
	"chr3R",
	"chr4",
	"chrX",
	"chrM",
	"chrU"
	]


GENES, EXONS, CDS, TRANSCRIPTS = {}, {}, {}, {}
exon_gene_ids, tx_gene_ids = set(), {}
gene_types = {}

out = open(args.output, "w")

for line in open(args.input):
	if line.startswith("##"):
		continue
	if line.strip() == "":
		continue
	line_el = line.strip().split("\t")
	if len(line_el) != 9:
		# The gff from FlyBase contains fasta lines as well
		continue
	chr, ann, element, start, end, score1, strand, score2, tags = line.strip().split("\t")
	if chr == "dmel_mitochondrion_genome":
		chr = "M"
	chr = "chr"+chr

	if ann != "FlyBase":
		continue
	if chr not in valid_chr:
		continue
	tags_d = dict(el.split("=") for el in tags.strip(";").split(";"))
	coords = "_".join((chr, start, end, strand))
	if element == "gene":
		GENES[tags_d["ID"]] = line
	if element == "exon":
		# Very few exons have the same name but different coordinates. 
		# Use both the coordinates and the name to have a unique identifier
		EXONS[";".join((coords, tags_d["Name"]))] = line
	if element == "CDS":
		CDS[";".join((coords, tags_d["Name"]))] = line
	if element in transcript_types:
		TRANSCRIPTS[tags_d["ID"]] = line
		gene_types.setdefault(tags_d["Parent"], []).append(element)

for exon,line in EXONS.iteritems():
	line_el = line.strip().split("\t")
	if line_el[0] == "dmel_mitochondrion_genome":
		line_el[0] = "M"
	line_el[0] = "chr" + line_el[0]
	tags_d = dict(el.split("=") for el in line_el[-1].strip(";").split(";"))
	transcripts = tags_d["Parent"].split(",")
	for transcript_id in transcripts:
		transcript_line = TRANSCRIPTS[transcript_id]
		transcript_el = transcript_line.strip().split("\t")
		tx_tags_d = dict(el.split("=") for el in transcript_el[-1].strip(";").split(";"))
		gene_id = tx_tags_d["Parent"]
		exon_gene_ids.add(gene_id)
		transcript_name = tx_tags_d["Name"]
		transcript_type = transcript_el[2]
		gene_type = transcript_types[max(transcript_types.index(x) for x in gene_types[gene_id])]
		gene_line = GENES[gene_id]
		gene_el = gene_line.strip().split("\t")
		gene_tags_d = dict(el.split("=") for el in gene_el[-1].strip(";").split(";"))
		gene_name = gene_tags_d["Name"]
		new_tag = 'gene_id "%s"; transcript_id "%s"; gene_name "%s"; transcript_name "%s"; gene_type "%s"; transcript_type "%s";' %(gene_id, transcript_id, gene_name, transcript_name, gene_type, transcript_type)
		out.write("\t".join(line_el[:-1]+ [new_tag]) + "\n")

for exon,line in CDS.iteritems():
	line_el = line.strip().split("\t")
	if line_el[0] == "dmel_mitochondrion_genome":
		line_el[0] = "M"
	line_el[0] = "chr" + line_el[0]
	tags_d = dict(el.split("=") for el in line_el[-1].strip(";").split(";"))
	transcripts = tags_d["Parent"].split(",")
	for transcript_id in transcripts:
		transcript_line = TRANSCRIPTS[transcript_id]
		transcript_el = transcript_line.strip().split("\t")
		tx_tags_d = dict(el.split("=") for el in transcript_el[-1].strip(";").split(";"))
		gene_id = tx_tags_d["Parent"]
		exon_gene_ids.add(gene_id)
		transcript_name = tx_tags_d["Name"]
		transcript_type = transcript_el[2]
		gene_type = transcript_types[max(transcript_types.index(x) for x in gene_types[gene_id])]
		gene_line = GENES[gene_id]
		gene_el = gene_line.strip().split("\t")
		gene_tags_d = dict(el.split("=") for el in gene_el[-1].strip(";").split(";"))
		gene_name = gene_tags_d["Name"]
		new_tag = 'gene_id "%s"; transcript_id "%s"; gene_name "%s"; transcript_name "%s"; gene_type "%s"; transcript_type "%s";' %(gene_id, transcript_id, gene_name, transcript_name, gene_type, transcript_type)
		out.write("\t".join(line_el[:-1]+ [new_tag]) + "\n")


for transcript, line in TRANSCRIPTS.iteritems():
	line_el = line.strip().split("\t")
	if line_el[0] == "dmel_mitochondrion_genome":
		line_el[0] = "M"
	line_el[0] = "chr" + line_el[0]
	transcript_type = line_el[2]
	line_el[2] = "transcript"
	tags_d = dict(el.split("=") for el in line_el[-1].strip(";").split(";"))
	transcript_name = tags_d["Name"]
	transcript_id = tags_d["ID"]
	gene_id = tags_d["Parent"]
	gene_line = GENES[gene_id]
	gene_el = gene_line.strip().split("\t")
	gene_tags_d = dict(el.split("=") for el in gene_el[-1].strip(";").split(";"))
	gene_name = gene_tags_d["Name"]
	gene_type = transcript_types[max(transcript_types.index(x) for x in gene_types[gene_id])]
	tx_gene_ids.setdefault(gene_id, []).append((transcript_id, transcript_name))
	new_tag = 'gene_id "%s"; transcript_id "%s"; gene_name "%s"; transcript_name "%s"; gene_type "%s"; transcript_type "%s";' %(gene_id, transcript_id, gene_name, transcript_name, gene_type, transcript_type)
	out.write("\t".join(line_el[:-1]+ [new_tag]) + "\n")


for gene, line in GENES.iteritems():
	line_el = line.strip().split("\t")
	if line_el[0] == "dmel_mitochondrion_genome":
		line_el[0] = "M"
	line_el[0] = "chr" + line_el[0]
	line_el[2] = "gene"
	tags_d = dict(el.split("=") for el in line_el[-1].strip(";").split(";"))
	gene_id = tags_d["ID"]
	gene_name = tags_d["Name"]
	try:
		gene_type = transcript_types[max(transcript_types.index(x) for x in gene_types.get(gene_id, None))]
	except:
		gene_type = "NA"
	transcript_id = gene_id
	transcript_name = gene_name
	transcript_type = gene_type
	new_tag = 'gene_id "%s"; transcript_id "%s"; gene_name "%s"; transcript_name "%s"; gene_type "%s"; transcript_type "%s";' %(gene_id, transcript_id, gene_name, transcript_name, gene_type, transcript_type)
	out.write("\t".join(line_el[:-1]+ [new_tag]) + "\n")
	if gene_id not in exon_gene_ids:
		if gene_id in tx_gene_ids:
			# Ideally I should insert a line for each transcript but now I take the first
			transcript_id, transcript_name = tx_gene_ids[gene_id][0]
			new_tag = 'gene_id "%s"; transcript_id "%s"; gene_name "%s"; transcript_name "%s"; gene_type "%s"; transcript_type "%s";' %(gene_id, transcript_id, gene_name, transcript_name, gene_type, transcript_type)
		line_el[2] = "exon"
		out.write("\t".join(line_el[:-1]+ [new_tag]) + "\n")
	if gene_id not in tx_gene_ids:
		line_el[2] = "transcript"
		out.write("\t".join(line_el[:-1]+ [new_tag]) + "\n")

out.close()

#print len(GENES)
#print len(exon_gene_ids)
#print len(tx_gene_ids)
#print set(GENES.keys()).difference(tx_gene_ids)
