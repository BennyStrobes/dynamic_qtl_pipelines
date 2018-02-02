args = commandArgs(trailingOnly=TRUE)
library(qvalue)


null_file = args[1] # Input file
qvalue_file = args[2] # output file

# Load in data
data <- read.table(null_file,header=TRUE)

qobj <- qvalue(p=data$pvalue)

qvalz <- qobj$qvalues


write.table(qvalz, file = qvalue_file, quote = FALSE, sep = "\t")