# script for running ccrepe on an OTU summary table from QIIME results

# load ccrepe library
library(ccrepe)

# set input file and output file
file = '/Users/jimbo/qiime/study_1001/taxonomy_summaries/study_1001_closed_reference_otu_table_L6.txt'
output = '/Users/jimbo/Desktop/working/study_1001_data/ccrepe_results_lvl6'

# read in a datafile of sample/taxa data
otu_abundance <- read.delim(file, header=TRUE)

# Reset row names and transpose matrix to form that CCREPE wants it in
otu_data <- otu_abundance[,-1]
rownames(otu_data) <- otu_abundance[,1]
otu_data.t <- t(otu_data)

# run ccrepe
data.ccrepe <- ccrepe(otu_data.t, min.subj = 10, verbose=TRUE, make.output.table = TRUE)

#Export data:   	
write.table(data.ccrepe$output.table, output)
