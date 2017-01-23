#!/usr/bin/env Rscript

##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library("optparse"))

###########
# VARIABLES
###########

options(stringsAsFactors=FALSE)


##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrices you want to analyze"),
make_option(c("-l", "--log10"), action="store_true", default=FALSE, help="apply the log10 [default=%default]"),
make_option(c("-p", "--pseudocount"), type="double", help="specify a pseudocount for the log [default=%default]", default=1e-04),
#make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--orth"), help="orthologous matching, no header."),
make_option(c("-k", "--keep_NA"), action="store_true", help="keep NAs, instead of replacing them by zero. [default=%default]", default=FALSE),
make_option(c("-O", "--output"), help="output file name [default=%default]", default='out.tsv')
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)

##-----------------------##
## BEGIN                 ##
##-----------------------##

# read inputs
matrices = strsplit(opt$input_matrix, ",")[[1]] 
orth = read.table(opt$orth, h=F, col.names=c('hs','mm'))

# merge the matrices over the orthologous genes
for ( i in 1:length(matrices) ) {
	m = read.table(matrices[i], h=T)
	if (grepl('human', basename(matrices[i]))) {merge_col = 'hs'}
	if (grepl('mouse', basename(matrices[i]))) {merge_col = 'mm'}
	orth_m = merge(orth, m, by.x=merge_col, by.y='row.names')
#	orth_m = merge(orth, m, by.x=merge_col, by.y='row.names', all.x=T)
	#print(matrices[i])
	#print(head(orth_m[,1:5]))
	if (i==1) {new_m = orth_m}
	else {new_m = merge(new_m, orth_m, by=c('hs','mm'))}
	#print(dim(new_m))
}

#print(head(new_m[,90:95]))

# modify the matrix according to user options
# -- replace NAs with 0
if (!opt$keep_NA) {new_m <- replace(new_m,is.na(new_m),0)}
# -- calculate the log and add the pseudocount
if (opt$log10) {new_m <- log10(new_m+opt$pseudocount)}

# compose the output file name
#output = sprintf("log10_%s.pdcn_%s.NAs_%s.%s.tsv", opt$log, opt$pseudocount, opt$keep_NA, opt$output)

output = opt$output
write.table(new_m, file=output, quote=F, sep="\t", row.names=F)

q(save='no')



