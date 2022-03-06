library("optparse")

# Rscript ../selectSamples.R --fam_file=${fam_file} --number=10000 --outFile1="selectedFam.fam"
# Output: 
# selectedFam.fam

option_list = list(
  make_option(c("-f", "--fam_file"), type = "character"),
  make_option(c("-n", "--number"), type = "numeric"),
  make_option(c("-s", "--outFile1"), type = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #==================================================
# opt$fam_file="qced_rdm_005_ukbb3.fam"
# opt$number=10000
# opt$outFile1="selectedFam.fam"
# #==================================================

#-------------------------------
# Select samples

fam <- read.table(opt$fam_file, stringsAsFactors=F, header=F)
total <- nrow(fam)

set.seed(2)
id <- sample(1:total, opt$number)
idv = fam[id,2]
idv <- sort(idv)

filename = opt$outFile1
writeLines(paste("Writing file: ", filename, sep = ""))
sink(filename)
write.table(fam[which(fam$V2 %in% idv),], row.names = F, quote = F, col.names = F)
sink()
#-------------------------------