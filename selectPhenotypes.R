library("optparse")

#Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --other_pheno_file=${other_pheno_file} --covariate_file=${covariate_file} --selected_id=${selected_id} --out_file=${out_file} --binary=${binary} --pre=${pre}

# Out
# out_file

option_list = list(
  make_option(c("--pheno_file"), type = "character"),
  make_option(c("--other_pheno_file"), type = "character"),
  make_option(c("--covariate_file"), type = "character"),
  make_option(c("--selected_id"), type = "character"),
  make_option(c("--out_file"), type = "character"),
  make_option(c("--binary"), type = "logical", default = F),
  make_option(c("--pre"), type = "character", default=""),
  make_option(c("--thres"), type = "numeric", default=20),
  make_option(c("--useSex"), type = "logical", default = T)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #==================================================
# opt$pheno_file="sh_trait.txt"
# opt$other_pheno_file="ukb6247.csv"
# opt$covariate_file="my_sample_qcinfo.txt"
# opt$selected_id="selectedFam.fam"
# opt$out_file="sh_random_50K_adjusted.dat"
# opt$pre="sh"
# #==================================================

#-------------------------------
# Read data
writeLines("Reading data ...")

fam=read.table(opt$selected_id)
head(fam)
nrow(fam)
# > head(fam)
# V1      V2 V3 V4 V5 V6
# 1 1000105 1000105  0  0  1 -9
# 2 1000137 1000137  0  0  2 -9
# 3 1000255 1000255  0  0  2 -9
# 4 1000323 1000323  0  0  2 -9
# 5 1000357 1000357  0  0  1 -9
# 6 1000439 1000439  0  0  1 -9
# > nrow(fam)
# [1] 50000

phe=read.table(opt$pheno_file)
colnames(phe) <- c("id", "trait")
head(phe)
nrow(phe)
# > head(phe)
# id edu
# 1 2080762  10
# 2 3165676  20
# 3 1605581  10
# 4 5336523  20
# 5 2211334  15
# 6 5391704  19
# > nrow(phe)
# [1] 502647

cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,$4,$3,$2353,$170,$18,$2354}' ", opt$other_pheno_file, sep = "") # eid, year of birth, sex, age, Townsend deprivation index at recruitment (TDI), centre, batch
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
# > head(dat)
# V1     V2     V3        V4        V5     V6        V7
# 1     eid 34-0.0 31-0.0 21022-0.0   189-0.0 54-0.0 22000-0.0
# 2 2080762   1958      1        49  -3.62553  11009      2000
# 3 3165676   1952      1        57 -0.402275  11020        15
# 4 1605581   1946      0        63  -4.00346  11017      2000
# 5 5336523   1953      1        54  -5.38678  11008      2000
# 6 2211334   1957      1        52 -0.248725  11009         3

cmd = paste("awk '{print $1,$2,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41}' ", opt$covariate_file, sep = "") # id & principal components
dat2 = read.table(pipe(cmd), header=T)
# > head(dat2)
# fid     iid      pc1     pc2       pc3       pc4      pc5       pc6
# 1 5438743 5438743 -12.1725 5.39163 -1.281030  0.841765 -5.26521 -1.786570
# 2 2890570 2890570 -13.0245 6.41514 -0.183365  2.927610 -5.88964  0.940534
# 3 2689958 2689958 -11.4712 3.48383 -1.154580  3.083830  7.65160 -0.913399
# 4 5800351 5800351 -12.1327 4.02976 -0.988080  0.750294 -2.36431  0.431658
# 5 5644632 5644632 -12.2171 3.50821 -1.625990 -1.226800 -5.34580  3.816790
# 6 5121015 5121015 -11.5847 5.41131  0.785815 -0.940348 -6.08543  1.465410
# pc7        pc8      pc9      pc10
# 1  3.109920 -2.6308500  2.39288  0.307537
# 2  1.141060 -1.9821300 -2.70226  2.507750
# 3 -1.548790  1.4789300 -1.20895  1.064900
# 4 -0.534071 -0.6543670 -6.59593 -1.533560
# 5  0.579155 -1.1695800  1.01855 -0.795227
# 6  2.698600  0.0793596 -2.85344  3.514010

#-------------------------------

#-------------------------------
# Check with fam file

# trait
m1 <- match(fam$V1, phe$id)
dat0 <- phe[m1,]
dat0 <- cbind(dat0[,1], dat0)
# > head(dat0)
# dat0[, 1]      id edu
# 460022   1000105 1000105  19
# 51796    1000137 1000137  20
# 245527   1000255 1000255  20
# 377680   1000323 1000323   7
# 241319   1000357 1000357  19
# 205486   1000439 1000439  19
# > nrow(dat0)
# [1] 50000

if(opt$binary) {
  dat0[,3] <- as.numeric(as.character(dat0[,3]))
  for (i in 1:nrow(dat0)) {
    if(!is.na(dat0[i,3])) {
      if (dat0[i,3] >= opt$thres) {
        dat0[i,3] <- 1
      } else {
        dat0[i,3] <- 0
      }
    }
  }
}
# > head(dat0)
# dat0[, 1]      id edu
# 460022   1000105 1000105   0
# 51796    1000137 1000137   1
# 245527   1000255 1000255   1
# 377680   1000323 1000323   0
# 241319   1000357 1000357   0
# 205486   1000439 1000439   0

# covariates
m1 <- match(fam$V2, dat$V1)
m2 <- match(fam$V2, dat2$iid)

year_of_birth <- as.numeric(as.character(dat[m1,2]))
sex <- as.character(dat[m1,3])
age <- as.numeric(as.character(dat[m1,4]))
tdi <- as.numeric(as.character(dat[m1,5]))
centre <- as.character(dat[m1,6])
batch <- as.character(dat[m1,7])
# > head(year_of_birth)
# [1] 1952 1964 1958 1951 1944 1949
# > length(year_of_birth)
# [1] 50000
# > head(sex)
# [1] "1" "0" "0" "0" "1" "1"
# > head(age)
# [1] 56 43 51 56 64 59
# > head(tdi)
# [1]  2.92534 -2.39857 -3.40437 -3.53948 -2.21228 -5.21046
# > head(centre)
# [1] "11009" "11003" "11021" "11007" "11010" "11013"
# > head(batch)
# [1] "2000" "10"   "2000" "2000" "2000" "2000"

sex <- paste("sex", sex, sep = "")
centre <- paste("centre", centre, sep = "")
batch <- paste("batch", batch, sep = "")
# > head(sex)
# [1] "sex1" "sex0" "sex0" "sex0" "sex1" "sex1"
# > head(centre)
# [1] "centre11009" "centre11003" "centre11021" "centre11007" "centre11010"
# [6] "centre11013"
# > head(batch)
# [1] "batch2000" "batch10"   "batch2000" "batch2000" "batch2000" "batch2000"

trait <- as.numeric(as.character(dat0[,3]))
# > head(trait)
# [1] 19 20 20  7 19 19
# > length(trait)
# [1] 50000

pc1 <- dat2[m2,3]
pc2 <- dat2[m2,4]
pc3 <- dat2[m2,5]
pc4 <- dat2[m2,6]
pc5 <- dat2[m2,7]
pc6 <- dat2[m2,8]
pc7 <- dat2[m2,9]
pc8 <- dat2[m2,10]
pc9 <- dat2[m2,11]
pc10 <- dat2[m2,12]
# > head(pc1)
# [1] -13.6236 -10.9322 -11.5181 -12.7950  -9.8953 -11.8684
# > length(pc1)
# [1] 50000

#-------------------------------

#-------------------------------
# Adjusting phenotypes

covar = data.frame(fam[, c(1,2)], trait, year_of_birth, sex, age, tdi, centre, batch, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)
rownames(covar) = fam$V2
# > head(covar)
# V1      V2 trait year_of_birth  sex age      tdi      centre
# 1000105 1000105 1000105    19          1952 sex1  56  2.92534 centre11009
# 1000137 1000137 1000137    20          1964 sex0  43 -2.39857 centre11003
# 1000255 1000255 1000255    20          1958 sex0  51 -3.40437 centre11021
# 1000323 1000323 1000323     7          1951 sex0  56 -3.53948 centre11007
# 1000357 1000357 1000357    19          1944 sex1  64 -2.21228 centre11010
# 1000439 1000439 1000439    19          1949 sex1  59 -5.21046 centre11013
# batch      pc1     pc2       pc3       pc4       pc5       pc6
# 1000105 batch2000 -13.6236 3.21643 -0.856048  0.750187 -2.017980 -0.350832
# 1000137   batch10 -10.9322 4.29625 -2.159000  0.156117 -6.480520  2.445700
# 1000255 batch2000 -11.5181 2.45262 -1.655760  0.412146  0.425523 -0.655667
# 1000323 batch2000 -12.7950 2.91524 -2.727940  3.615550  3.929570 -2.938990
# 1000357 batch2000  -9.8953 4.26850 -1.368900  5.141510 -5.555260 -1.006820
# 1000439 batch2000 -11.8684 4.13917 -1.176380 -0.254564 -9.686770  0.671343
# pc7      pc8       pc9      pc10
# 1000105  5.101410  2.18070 -6.043430  1.789280
# 1000137  1.729390 -2.68435 -2.295680  1.192040
# 1000255 -2.218970  1.54612  6.145580 -0.896603
# 1000323 -0.454737  2.31013  2.517830 -4.155920
# 1000357 -0.533635 -1.79225 -4.053090  0.770221
# 1000439 -2.784570  1.22552 -0.603017 -2.051650
# > nrow(covar)
# [1] 50000
#covar = subset(covar, complete.cases(covar))
colnames(covar) = c("FID", "IID", "trait", "year", "sex", "age", "tdi", "centre", "batch","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")

if(opt$useSex == FALSE) {
  covar <- covar[,-5]
}

# Adjusting phenotypes
if(opt$binary) {
  #model = glm(as.formula(paste("trait", "~.", sep ="")), data=covar[, 3:ncol(covar)], family = binomial)
  model = glm(as.formula(paste("trait", "~.", sep ="")), data=covar[, 3:ncol(covar)], family = gaussian)
} else {
  model = glm(as.formula(paste("trait", "~.", sep ="")), data=covar[, 3:ncol(covar)], family = gaussian)  
}
trait_adjusted = data.frame(names(model$residuals), names(model$residuals), model$residuals)
colnames(trait_adjusted) = c("FID", "IID", "trait_adjusted")
# > head(trait_adjusted)
# FID     IID trait_adjusted
# 1000105 1000105 1000105       4.820639
# 1000137 1000137 1000137       3.548171
# 1000255 1000255 1000255       4.579520
# 1000323 1000323 1000323      -7.977612
# 1000357 1000357 1000357       4.713562
# 1000439 1000439 1000439       3.141979
# > nrow(trait_adjusted)
# [1] 49471

trait_scale = data.frame(names(model$residuals), names(model$residuals), scale(model$residuals))
colnames(trait_scale) = c("FID", "IID", "trait_scale")
# > head(trait_scale)
# FID     IID trait_scale
# 1000105 1000105 1000105   0.9993857
# 1000137 1000137 1000137   0.7355853
# 1000255 1000255 1000255   0.9493983
# 1000323 1000323 1000323  -1.6538702
# 1000357 1000357 1000357   0.9771872
# 1000439 1000439 1000439   0.6513761
# > nrow(trait_scale)
# [1] 49471

#-------------------------------

#-------------------------------
# Writing files

writeLines("Writing files ...")

pre <- opt$pre
if(pre != "") {
  pre <- paste(pre, "_", sep = "")
}

if(opt$binary) {
  filename = paste(pre,"covariate_bi.tsv", sep ="")
} else {
  filename = paste(pre,"covariate.tsv", sep ="")
}
writeLines(paste("Writing file: ", filename, sep = ""))
sink(filename)
write.table(covar,quote=F,col.name=T,row.name=F, sep = "\t")
sink()

writeLines(paste("Writing file: ", opt$out_file, sep = ""))
sink(opt$out_file)
write.table(trait_adjusted,quote=F,col.name=F,row.name=F, sep = "\t")
sink()

# writeLines(paste("Writing file: ", "qualification_random_50K_adjusted_scale.dat", sep = ""))
# sink("qualification_random_50K_adjusted_scale.dat")
# write.table(trait_scale,quote=F,col.name=F,row.name=F, sep = "\t")
# sink()
#-------------------------------

res=model$residuals
head(res)

#covar = data.frame(fam[, c(1,2)], trait, year_of_birth, sex, age, tdi, centre, batch, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)
sv1=seq(1,length(fam$V1)) #50000
#sv2=sv1[!is.na(trait) & !is.na(sex) & !is.na(BY) & !is.na(TDI) & !is.na(centre) & !is.na(EDU) & !is.na(geno_batch)]
sv2=sv1[!is.na(trait) & !is.na(sex) & !is.na(year_of_birth) & !is.na(tdi) & !is.na(centre) & !is.na(age) & !is.na(batch)] #49471
v1=res #49471
out=array(NA,length(fam$V1)) #50000
out[sv2]=v1
if(opt$binary) {
  sink(paste(pre,"trait_adjusted_bi_f",sep=""))
} else {
  sink(paste(pre,"trait_adjusted_f",sep=""))  
}
write.table(cbind(as.character(fam$V1),as.character(fam$V2),out), quote=F, col.name=F, row.name=F)
sink()
