
  library(MASS)

  dat=read.table("selectedFam.fam")
  yv=read.table("mvnorm.out")

  out=yv$V1+yv$V2

  v1=cbind(as.character(dat$V1),as.character(dat$V2),out,yv$V1,yv$V2)
  sink("sample.dat")
  write.table (v1,row.names=F,col.names=F,quote=F)
  sink()


