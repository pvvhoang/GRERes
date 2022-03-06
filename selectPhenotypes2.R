# From Momin
#-------------------------------
setwd("C:/Users/phamvv/MyDoc/GRERes/Data")


getwd()

#Fam file ethnicity ##Final Fam file



#fam=read.table("C:/Users/mommy003/Desktop/Fam file ethnicity/ethnicity1b.fam", header = F)
fam=read.table("selectedFam.fam", header = F)

dim(fam)

head(fam)





#covariate=read.table("covariate",header=T) 
covariate=read.table("dat_covariate",header=T) 

head(covariate)

colnames(covariate)

#colnames(covariate)=c("eid","age","sex","BY","TDI","centre","geno_batch")
colnames(covariate)=c("eid","BY","sex","age","TDI","centre","geno_batch")

head(covariate)

dim(covariate)







mat=match(fam$V1,covariate$eid)

length(mat)

covar=covariate[mat,]

dim(covar) ##final number for adjustment

head(covar)





#pc=read.table("pc_1_10", header = T)
pc=read.table("dat2_pc", header = T)

head(pc)

dim(pc)





# id alignment

m=match(covar$eid,pc$fid)

dat=data.frame(id=covar$eid,
               
               pc[m,-c(1,2)])

length (dat)

head(dat)

dim(dat)





# trait=read.table("2. standing_height",header=T)
# 
# head(trait)
# 
# summary(trait)
# 
# dim(trait)
# 
# m2=match(fam$V1,trait$eid)
# 
# length(m2)
# 
# trait=trait[m2,]
# 
# dim(trait)





edu_qua=read.table("edu_corrected")

head(edu_qua)

m1=match(fam$V1,edu_qua$id)

EDU=edu_qua[m1,]

head(EDU)

dim(EDU)





#trait=trait$X50.0.0
trait=EDU$edu

summary(trait)

age=covar$age

sex=covar$sex

BY=covar$BY

TDI=covar$TDI

#EDU=EDU$edu

centre=as.factor(covar$centre)

geno_batch=as.factor(covar$geno_batch)







pc1=dat$pc1

pc2=dat$pc2

pc3=dat$pc3

pc4=dat$pc4

pc5=dat$pc5

pc6=dat$pc6

pc7=dat$pc7

pc8=dat$pc8

pc9=dat$pc9

pc10=dat$pc10









#adj=lm(trait~age+sex+BY+TDI+centre+EDU+geno_batch+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10)
adj=lm(trait~age+sex+BY+TDI+centre+geno_batch+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10)





summary(adj)

# > summary(adj)
# 
# Call:
#   lm(formula = trait ~ age + sex + BY + TDI + centre + geno_batch + 
#        pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -12.396  -4.702   1.278   4.264   9.999 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    -3.271e+02  8.565e+01  -3.819 0.000134 ***
#   age             3.267e-02  4.269e-02   0.765 0.444110    
# sex             9.605e-01  4.362e-02  22.023  < 2e-16 ***
#   BY              1.738e-01  4.268e-02   4.073 4.65e-05 ***
#   TDI            -2.722e-01  8.017e-03 -33.954  < 2e-16 ***
#   centre11002     7.659e-01  1.820e-01   4.208 2.58e-05 ***
#   centre11003    -4.522e-01  1.807e-01  -2.503 0.012333 *  
#   centre11004     2.379e-01  1.782e-01   1.335 0.181798    
# centre11005     1.162e+00  1.784e-01   6.510 7.58e-11 ***
#   centre11006    -9.175e-01  1.767e-01  -5.194 2.07e-07 ***
#   centre11007     4.864e-02  1.667e-01   0.292 0.770404    
# centre11008    -8.085e-01  1.629e-01  -4.963 6.95e-07 ***
#   centre11009    -5.400e-01  1.624e-01  -3.325 0.000885 ***
#   centre11010    -5.172e-01  1.594e-01  -3.245 0.001177 ** 
#   centre11011    -5.960e-02  1.662e-01  -0.359 0.719868    
# centre11012     3.196e+00  2.255e-01  14.170  < 2e-16 ***
#   centre11013    -5.790e-01  1.686e-01  -3.435 0.000593 ***
#   centre11014    -3.393e-01  1.904e-01  -1.782 0.074767 .  
# centre11016    -4.555e-01  1.797e-01  -2.535 0.011251 *  
#   centre11017    -5.018e-01  1.913e-01  -2.623 0.008731 ** 
#   centre11018     9.703e-01  1.952e-01   4.970 6.71e-07 ***
#   centre11020     4.098e-01  1.989e-01   2.060 0.039371 *  
#   centre11021    -7.069e-01  2.002e-01  -3.531 0.000414 ***
#   centre11022    -6.472e-01  3.473e-01  -1.864 0.062378 .  
# centre11023    -1.419e+00  6.122e-01  -2.318 0.020441 *  
#   geno_batch-10  -2.398e-01  3.125e-01  -0.767 0.442926    
# geno_batch-9   -2.322e-01  3.160e-01  -0.735 0.462497    
# geno_batch-8   -3.708e-01  3.058e-01  -1.212 0.225360    
# geno_batch-7   -7.294e-02  3.059e-01  -0.238 0.811554    
# geno_batch-6   -2.252e-01  3.093e-01  -0.728 0.466617    
# geno_batch-5   -2.117e-01  3.028e-01  -0.699 0.484369    
# geno_batch-4    2.441e-02  3.079e-01   0.079 0.936812    
# geno_batch-3    1.948e-01  3.079e-01   0.633 0.526865    
# geno_batch-2   -5.858e-01  3.082e-01  -1.901 0.057344 .  
# geno_batch-1   -1.861e-01  3.084e-01  -0.604 0.546142    
# geno_batch1     3.565e-01  3.077e-01   1.159 0.246592    
# geno_batch2     1.112e-01  3.265e-01   0.340 0.733482    
# geno_batch3    -4.739e-02  3.167e-01  -0.150 0.881058    
# geno_batch4     5.105e-01  3.114e-01   1.640 0.101072    
# geno_batch5     4.745e-01  3.060e-01   1.551 0.120977    
# geno_batch6     1.416e-01  3.114e-01   0.455 0.649199    
# geno_batch7     4.490e-01  3.110e-01   1.444 0.148824    
# geno_batch8     3.784e-01  3.111e-01   1.216 0.223832    
# geno_batch9     2.757e-01  3.128e-01   0.881 0.378102    
# geno_batch10    3.382e-01  3.150e-01   1.074 0.282985    
# geno_batch11    2.654e-01  3.081e-01   0.862 0.388902    
# geno_batch12    1.803e-01  3.122e-01   0.578 0.563564    
# geno_batch13    1.186e-01  3.195e-01   0.371 0.710408    
# geno_batch14    5.230e-01  3.196e-01   1.636 0.101750    
# geno_batch15    1.382e-02  3.089e-01   0.045 0.964311    
# geno_batch16    4.439e-01  3.138e-01   1.415 0.157177    
# geno_batch17    4.216e-01  3.137e-01   1.344 0.179038    
# geno_batch18    2.217e-01  3.141e-01   0.706 0.480460    
# geno_batch19    4.027e-01  3.188e-01   1.263 0.206540    
# geno_batch20    4.317e-01  3.126e-01   1.381 0.167274    
# geno_batch21    1.131e+00  3.168e-01   3.569 0.000359 ***
#   geno_batch22    7.711e-02  3.185e-01   0.242 0.808711    
# geno_batch2000  2.709e-01  2.211e-01   1.225 0.220419    
# pc1            -7.897e-03  1.424e-02  -0.555 0.579208    
# pc2            -1.413e-03  1.474e-02  -0.096 0.923654    
# pc3            -1.315e-03  1.428e-02  -0.092 0.926596    
# pc4            -4.358e-02  1.064e-02  -4.097 4.20e-05 ***
#   pc5            -2.224e-03  4.898e-03  -0.454 0.649710    
# pc6            -5.656e-03  1.355e-02  -0.417 0.676494    
# pc7            -1.572e-02  1.212e-02  -1.297 0.194524    
# pc8            -3.043e-03  1.209e-02  -0.252 0.801205    
# pc9            -2.140e-02  5.520e-03  -3.876 0.000106 ***
#   pc10            3.238e-02  1.102e-02   2.939 0.003293 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 4.827 on 49403 degrees of freedom
# (529 observations deleted due to missingness)
# Multiple R-squared:  0.09064,	Adjusted R-squared:  0.0894 
# F-statistic: 73.49 on 67 and 49403 DF,  p-value: < 2.2e-16



res=adj$residuals

head(res)



sv1=seq(1,length(fam$V1)) #50000

#sv2=sv1[!is.na(trait) & !is.na(sex) & !is.na(BY) & !is.na(TDI) & !is.na(centre) & !is.na(EDU) & !is.na(geno_batch)]
sv2=sv1[!is.na(trait) & !is.na(sex) & !is.na(BY) & !is.na(TDI) & !is.na(centre) & !is.na(age) & !is.na(geno_batch)] #49471



v1=scale(adj$residuals) #49471

out=array(NA,length(fam$V1)) #50000

out[sv2]=v1

sink("trait_adjusted")

write.table(cbind(as.character(fam$V1),as.character(fam$V2),out), quote=F, col.name=F, row.name=F)

sink()



tmp <- cbind(as.character(fam$V1),as.character(fam$V2),out)
tmp = subset(tmp, complete.cases(tmp))

sink("trait_adjusted2")
write.table(tmp, quote=F, col.name=F, row.name=F)
sink()
