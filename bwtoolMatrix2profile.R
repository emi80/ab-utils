#!/usr/bin/env Rscript


options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c('-i', '--input'), type="character", default="stdin",
	help='Input matrix, tab-separated. [default=%default]'),

#make_option(c('-r', '--row_names'), type="numeric", help='Index of the column with the row names. Leave empty for no row names.'),
make_option(c('-o', '--output'), type="character", default='stdout',
	help='Output file name with extension. Can be <stdout>. [default=%default]'),

make_option(c("--header"), action="store_true", default=FALSE,
	help="The input matrix has a header [default=%default]"),

make_option(c("-r", "--replace_NAs"), action="store_true", default=FALSE,
	help="Replace NAs by 0s [default=%default]"),

make_option(c('-a', '--aggregate_by'), type="numeric", 
	help='Index of the column with the factor to aggregate by. Leave empty for no factor.'),

make_option(c('-f', '--func'), type="character", default='mean', 
	help='Aggregating function [default=%default]'),

make_option(c("-B", "--no_reshape"), action="store_true", default=FALSE,
	help="Do not reshape the data frame [default=%default]"),

make_option(c('-n', '--offset'), type="numeric", default=0, help='Offset [default=%default]'),
make_option(c('-v', '--verbose'), action='store_true', default=FALSE, help='Verbose output [default=%default]')

)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


#~~~~~~~~~~~~~~~~~~~~~
# LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~



if (opt$verbose) {cat("Loading libraries... ")}
suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}





#########
# BEGIN
#########

# read input
if(opt$input == 'stdin') {input=file('stdin')} else {input=file(opt$input)}
m = read.table(input, h=opt$header)

# Replace NAs
if (opt$replace_NAs) {m=replace(m, is.na(m), 0)}

# Read function
func = eval(opt$func)

# Aggregate and apply function
if (!is.null(opt$aggregate_by)) {
	# read factor index
	factor_col = colnames(m)[opt$aggregate_by]
	form = as.formula(sprintf(".~%s", factor_col))
	agg = aggregate(form, m, func, na.rm=T)
} else {
	agg = apply(m, 2, func, na.rm=T)
}


# Reshape if needed

if (!opt$no_reshape) {
	# Melt the data
	if (!is.null(opt$aggregate_by)) {
		mm = melt(agg, id.vars=factor_col, variable.name=c('position'))
		# Reorder column
		mm <- mm[c(2,3,1)]
	} else {
		mm = data.frame(position=names(agg), value=agg)
	}
	
	# Read the offset and edit the position column accordingly
	offset = opt$offset
	new_lev = seq_along(levels(mm$position)) - offset
	levels(mm$position) <- new_lev
} else {
	mm = agg
}

#~~~~~~~~~~~~
# OUTPUT
#~~~~~~~~~~~~

output = ifelse(opt$output=='stdout', "", opt$output)
write.table(mm, file=output, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

q(save='no')


# !! TO DO: !!
# - Read more than one column as a factor
# - Compute more than one function on the data

