library("optparse")
library(MASS)

# Rscript ../toBinary.R --phe="sample_pheno.dat"
# Output: 
# sample_pheno.dat

option_list = list(
  make_option(c("-p", "--phe"), type = "character"),
  make_option(c("-k", "--kval"), type = "numeric", default = 0.1)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #=================================================
# opt$phe="sample_pheno.dat"
# #=================================================
  
  # case-control data with k=0.1
  # dnorm (thd) = height at the the threshold in the standard normal curve.
  # -qnorm (k) = thd
  # pnorm (thd) = k
  k <- opt$kval
  thd <- -qnorm (k)
  
  phe <- read.table(opt$phe, stringsAsFactors=F, header=F)
  # > head(phe)
  # V1      V2          V3
  # 1 1000584 1000584  0.05143439
  # 2 1001316 1001316  0.58217016
  # 3 1001492 1001492  0.68453999
  # 4 1002199 1002199  2.69808509
  # 5 1002524 1002524 -0.26664219
  # 6 1002548 1002548  0.08108904
  # > mean(phe[,3])
  # [1] -3.005185e-17
  # > sd(phe[,3])
  # [1] 1.003653
  phe[,3] <- scale(phe[,3])
  # > head(phe)
  # V1      V2          V3
  # 1 1000584 1000584  0.05124720
  # 2 1001316 1001316  0.58005144
  # 3 1001492 1001492  0.68204870
  # 4 1002199 1002199  2.68826578
  # 5 1002524 1002524 -0.26567178
  # 6 1002548 1002548  0.08079393
  # > mean(phe[,3])
  # [1] -1.490346e-17
  # > sd(phe[,3])
  # [1] 1
  phe[,4] <- 0
  idx <- which(phe[,3] > thd)
  phe[idx,4] <- 1
  phe <- phe[,-3]
  # >   ind <- which(phe[,3] == 1)
  # >   length(ind)
  # [1] 972
  
  sink("sample_pheno.dat")
  write.table(phe,quote=F,col.name=F,row.name=F)
  sink()
  