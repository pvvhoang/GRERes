library(sm)

isBiased <- function(estimatedValues, trueValue) {
  sel <- estimatedValues
  ave <- mean(sel, na.rm=T)
  std <- sd(sel, na.rm=T)
  if ((trueValue < ave+(1.96*std)) & (trueValue > ave-(1.96 * std))) {
    return(FALSE)
  }
  return(TRUE)
}

drawOneChart <- function(label="Covariance", dat, file) {
  # pdf(file = file,   # The directory you want to save the file in
  #     width = 6, # The width of the plot in inches
  #     height = 6) # The height of the plot in inches
  jpeg(file = file,   # The directory you want to save the file in
       width = 350, # The width of the plot in inches
       height = 350) # The height of the plot in inches
  
  # Filled Density Plot
  d <- density(dat)
  plot(d, main=label)
  polygon(d, col="#ADD8E6", border="#ADD8E6")
  
  dev.off()
  
}

drawTwoCharts <- function(label1="GREML", dat1, label2="COREGREML", dat2, file, title="") {
  # pdf(file = file,   # The directory you want to save the file in
  #     width = 6, # The width of the plot in inches
  #     height = 6) # The height of the plot in inches
  jpeg(file = file,   # The directory you want to save the file in
       width = 350, # The width of the plot in inches
       height = 350) # The height of the plot in inches
  
  n1 <- length(dat1)
  n2 <- length(dat2)
  d <- matrix(0, nrow=n1+n2, ncol=2)
  d[1:n1,1] <- 1 # for label1
  d[(n1+1):(n1+n2),1] <- 2 # for label2
  d[1:n1,2] <- dat1
  d[(n1+1):(n1+n2),2] <- dat2

  l.f <- factor(d[,1], levels= c(1, 2),                               
                  labels = c(label1, label2)) 
  
  sm.density.compare(d[,2], d[,1], xlab="Estimate")                
  title(main=title)
  
  colfill<-c(2:(2+length(levels(l.f)))) 
  
  legend("topleft",locator(1), levels(l.f), fill=colfill) 
  
  dev.off()
}

computeType1Error <- function(cov, se) {
  n <- length(cov)
  type1_rate <- c()
  for (i in 1:n) {
    chi = (cov[i]/se[i])^2
    p <- pchisq(chi,1,lower.tail=F)
    type1_rate <- c(type1_rate, p)
  }
  idx <- which(type1_rate < 0.05)
  # percentage of type 1 error rate
  per <- length(idx)/length(type1_rate)*100
  return(per)
}
