#!/usr/bin/env python
import os
import csv
import sys
import argparse
from xlsxwriter.workbook import Workbook

class ParseLabels(argparse.Action):
     def __init__(self, option_strings, dest, nargs=None, **kwargs):
         if nargs is not None:
             raise ValueError("nargs not allowed")
         super(ParseLabels, self).__init__(option_strings, dest, **kwargs)
     def __call__(self, parser, namespace, values, option_string=None):
         d = {}
         with open(values) as f:
         	for line in f:
         		l = line.strip("\n").split("\t")
         		d[l[0]] = l[1]
         setattr(namespace, self.dest, d)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Add tsv file(s) to excel')
	parser.add_argument('tsv_files', type=str, nargs='+',
                    help='input files')
	parser.add_argument('-o', '--out-file', dest='outfile', 
	                    help='the output file (with or wihout the extension)')
	parser.add_argument('-l', '--labels-file', dest='labels', action=ParseLabels, 
	                    help='TSV file containing the mapping between file name and sheet label')

	args = parser.parse_args()

	# Add some command-line logic to read the file names.
	xlsx_file = args.outfile

	if '.xlsx' not in xlsx_file:
		xlsx_file = "{}.xlsx".format(xlsx_file)

	# Create an XlsxWriter workbook object and add a worksheet.
	workbook = Workbook(xlsx_file)

	for i, tsv_file in enumerate(args.tsv_files):
		label = args.labels.get(os.path.basename(tsv_file),"sheet{}".format(i))
		worksheet = workbook.add_worksheet(label)

		# Create a TSV file reader.
		tsv_reader = csv.reader(open(tsv_file, 'rb'), delimiter='\t')

		# Read the row data from the TSV file and write it to the XLSX file.
		for row, data in enumerate(tsv_reader):
		    worksheet.write_row(row, 0, data)

	# Close the XLSX file.
	workbook.close()