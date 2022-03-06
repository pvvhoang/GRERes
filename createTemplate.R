library("optparse")
library(MASS)

# Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
# Output: 
# outFile

option_list = list(
  make_option(c("-i", "--inFile"), type = "character"),
  make_option(c("-o", "--outFile"), type = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #=================================================
# opt$inFile="brain_coregreml.out2"
# opt$outFile="brain_coregreml.do"
# #=================================================
  
sink(opt$outFile)
cat(opt$inFile, "\n") # line 1: specify the file with parameter estimates
cat(3, "\n") # line 2: tot. # of variance & covariance components in the file
cat("R 2 1 3 3", "\n") # line 3: compute prop. of variance due to genetics ('R' is to get the ratio of variance #2, i.e. #2 / (#1 + #2 + #3 + #3))
cat("R 1 2 3 3", "\n")
cat("C 3 1 2", "\n") # line 5: compute correlation between g & b ('C' is to get the correlation of covariance #3, i.e. #3 / sqrt(#1 * #2))
sink()
  