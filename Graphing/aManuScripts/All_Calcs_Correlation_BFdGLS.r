
rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library(xlsx)

##################################################################################################
## Test correlation across all hypotheses 
##################################################################################################
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# https://rpubs.com/aaronsc32/spearman-rank-correlation

emptyFrame <- function(h,s){
  c.d <- data.frame("test"=c("emp.spearman","emp.pearson"),
                          "sample"=c(s,s),
                          "hypothesis"=c(h,h),
                          "Correlation"=c(NA,NA),
                          "emp"=c(NA,NA),
                          "p0.975"=c(NA,NA),
                          "p0.95"=c(NA,NA),
                          "p0.05"=c(NA,NA),
                          "p0.025"=c(NA,NA))
  return(c.d)
}

corTest <- function(d,h,s){
  
  # Create two datasets, one to add permutation columns to, one to pull from. because that's how I coded it?
  a <- d %>% filter(variable == h) 
  perms <- d %>% filter(variable == h) 
  nloci <- length(a$Locus)
  maxB <- max(a$BF)
  minB <- min(a$BF)
  maxG <- max(a$GL)
  minG <- min(a$GL)
  #print(length(a$Locus))
  
  # Safety check for if a subset of data has too few observations
  if(length(a$Locus) < 5){
    c.d <- emptyFrame(h,s)
    return(c.d)
  }
  
  # Define which column numbers to test for correlation
  
  b <- match('BF',names(a))
  g <- match('GL',names(a))

  # Permute BF values in 1,000 seperate columns, could permute BF or GL doesnt matter. 
  
  for (n in 1:1000) {
    newcolname1 <- paste('BF',n,sep='') 
    perms[, newcolname1] <- sample(a[[b]], replace = FALSE)
    #newcolname2 <- paste('GL',n,sep='') 
    #null[, newcolname2] <- sample(a[[g]], replace = FALSE)
  }
  
  # Create new dataset to add null rank coefficients to 
  
  assign('null', setNames(data.frame(1:1000), "replicate"))
  
  # Add bf gl 
  
  for (n in 1:1000) {
    
    # Grab column names, calculate value 
    colname1 <- paste('BF',n,sep='') 
    cc1.s <- cor(perms[,'GL'],perms[,colname1],method="spearman") 
    null[n,'null.spearman'] <- cc1.s
    cc1.p <- cor(perms[,'GL'],perms[,colname1],method="pearson") 
    null[n,'null.pearson'] <- cc1.p
  }
  
  
  # Dataframe for emp scores
  emp <- data.frame(TestType=character(),
                    corr.eff=integer(),
                    p=integer(),
                    p0.975=integer(),
                    p0.95=integer(),
                    p0.05=integer(),
                    p0.025=integer(),stringsAsFactors=FALSE)
  
  # Designate which row is for which test 
  
  emp[1,'TestType'] <- 'emp.spearman'
  emp[2,'TestType'] <- 'emp.pearson'
  

  # Calculate correlation stats from emperical data,  add to df

  e.s <- cor(a[[b]],a[[g]],method="spearman") 
  e.p <- cor(a[[b]],a[[g]],method="pearson") 
  
  emp[1,'corr.eff'] <- e.s
  emp[2,'corr.eff'] <- e.p
  
  
  # Calculate 2 sided p values
  
  p.s <- sum(abs(null$null.spearman) >= abs(e.s)) / length(null$null.spearman)
  p.p <- sum(abs(null$null.pearson) >= abs(e.s)) / length(null$null.pearson)
  
  emp[1,'p'] <- p.s
  emp[2,'p'] <- p.p
  
  
  # Add rows for quantile cut offs 
  
  emp[1,'p0.025'] <- quantile(null$null.spearman, c(0.025, 0.975),na.rm = TRUE)[[1]]
  emp[1,'p0.975'] <- quantile(null$null.spearman, c(0.025, 0.975),na.rm = TRUE)[[2]]
  emp[1,'p0.05'] <- quantile(null$null.spearman, c(0.05, 0.95),na.rm = TRUE)[[1]]
  emp[1,'p0.95'] <- quantile(null$null.spearman, c(0.05, 0.95),na.rm = TRUE)[[2]]
  emp[2,'p0.025'] <- quantile(null$null.pearson, c(0.025, 0.975),na.rm = TRUE)[[1]]
  emp[2,'p0.975'] <- quantile(null$null.pearson, c(0.025, 0.975),na.rm = TRUE)[[2]]
  emp[2,'p0.05'] <- quantile(null$null.pearson, c(0.05, 0.95),na.rm = TRUE)[[1]]
  emp[2,'p0.95'] <- quantile(null$null.pearson, c(0.05, 0.95),na.rm = TRUE)[[2]]
  
  # Add significance yes/no and other data 
  
  corr.data <- emp %>%
    mutate_if(is.numeric,round,2) %>% 
    add_column(Hypothesis = h, .after = "TestType") %>% 
    add_column(Sample = s, .after = "TestType") %>% 
    add_column(Loci = nloci, .after = "p") %>%  
    add_column(minG = minG, .after = "Loci")%>% 
    add_column(maxG = maxG, .after = "Loci")%>% 
    add_column(minB = minB, .after = "Loci")%>% 
    add_column(maxB = maxB, .after = "Loci")
  
  return(corr.data)
}

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/Tables")



tvs <- corTest(all,"TvS","all") 
aivsa <- corTest(all,"AIvSA","all")
aivsi <- corTest(all,"AIvSI","all")
savsi <- corTest(all,"SAvSI","all")

tvs.0x <- corTest(all.0x,"TvS","sig0x")
aivsa.0x <- corTest(all.0x,"AIvSA","sig0x")
aivsi.0x <- corTest(all.0x,"AIvSI","sig0x")
savsi.0x <- corTest(all.0x,"SAvSI","sig0x")

tvs.1x <- corTest(all.1x,"TvS","sig1x")
aivsa.1x <- corTest(all.1x,"AIvSA","sig1x")
aivsi.1x <- corTest(all.1x,"AIvSI","sig1x")
savsi.1x <- corTest(all.1x,"SAvSI","sig1x")

tvs.2x <- corTest(all.2x,"TvS","sig2x")
aivsa.2x <- corTest(all.2x,"AIvSA","sig2x")
aivsi.2x <- corTest(all.2x,"AIvSI","sig2x")
savsi.2x <- corTest(all.2x,"SAvSI","sig2x")


tvs.50m <- corTest(all.50m,"TvS","quant50m")
aivsa.50m <- corTest(all.50m,"AIvSA","quant50m")
aivsi.50m <- corTest(all.50m,"AIvSI","quant50m")
savsi.50m <- corTest(all.50m,"SAvSI","quant50m")


tvs.60m <- corTest(all.60m,"TvS","quant60m")
aivsa.60m <- corTest(all.60m,"AIvSA","quant60m")
aivsi.60m <- corTest(all.60m,"AIvSI","quant60m")
savsi.60m <- corTest(all.60m,"SAvSI","quant60m")


tvs.80m <- corTest(all.80m,"TvS","quant80m")
aivsa.80m <- corTest(all.80m,"AIvSA","quant80m")
aivsi.80m <- corTest(all.80m,"AIvSI","quant80m")
savsi.80m <- corTest(all.80m,"SAvSI","quant80m")

tvs.60mb <- corTest(all.60mb,"TvS","quant60mb")
aivsa.60mb <- corTest(all.60mb,"AIvSA","quant60mb")
aivsi.60mb <- corTest(all.60mb,"AIvSI","quant60mb")
savsi.60mb <- corTest(all.60mb,"SAvSI","quant60mb")


tvs.80mb <- corTest(all.80mb,"TvS","quant80mb")
aivsa.80mb <- corTest(all.80mb,"AIvSA","quant80mb")
aivsi.80mb <- corTest(all.80mb,"AIvSI","quant80mb")
savsi.80mb <- corTest(all.80mb,"SAvSI","quant80mb")


tvs.50l <- corTest(all.50l,"TvS","50l")
aivsa.50l <- corTest(all.50l,"AIvSA","50l")
aivsi.50l <- corTest(all.50l,"AIvSI","50l")
savsi.50l <- corTest(all.50l,"SAvSI","50l")

tvs.50h <- corTest(all.50h,"TvS","50h")
aivsa.50h <- corTest(all.50h,"AIvSA","50h")
aivsi.50h <- corTest(all.50h,"AIvSI","50h")
savsi.50h <- corTest(all.50h,"SAvSI","50h")



results <- bind_rows(tvs,aivsa,aivsi,savsi,
                     tvs.80m,aivsa.80m,aivsi.80m,savsi.80m,
                     tvs.60m,aivsa.60m,aivsi.60m,savsi.60m,
                     tvs.50m,aivsa.50m,aivsi.50m,savsi.50m,
                     tvs.50l,aivsa.50l,aivsi.50l,savsi.50l,
                     tvs.50h,aivsa.50h,aivsi.50h,savsi.50h,
                     tvs.2x,aivsa.2x,aivsi.2x,savsi.2x,
                     tvs.1x,aivsa.1x,aivsi.1x,savsi.1x,
                     tvs.0x,aivsa.0x,aivsi.0x,savsi.0x) %>% filter(TestType == "emp.spearman")

results_mb <- bind_rows(tvs,aivsa,aivsi,savsi,
                     tvs.80m,aivsa.80m,aivsi.80m,savsi.80m,
                     tvs.60m,aivsa.60m,aivsi.60m,savsi.60m,
                     tvs.80mb,aivsa.80mb,aivsi.80mb,savsi.80mb,
                     tvs.60mb,aivsa.60mb,aivsi.60mb,savsi.60mb
                     ) %>% filter(TestType == "emp.spearman")


library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/Tables/AllSchmoosh.xlsx",sep='')

write.xlsx2(as.data.frame(results), file=fname, row.names=FALSE)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/Tables/AllSchmoosh_mb.xlsx",sep='')

write.xlsx2(as.data.frame(results_mb), file=fname, row.names=FALSE)


##################################################################################################
## Test correlation blind to hypotheses 
##################################################################################################
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# https://rpubs.com/aaronsc32/spearman-rank-correlation

emptyFrame <- function(s){
  c.d <- data.frame("test"=c("emp.spearman","emp.pearson"),
                    "sample"=c(s,s),
                    "hypothesis"=c('ALL','ALL'),
                    "Correlation"=c(NA,NA),
                    "emp"=c(NA,NA),
                    "p0.975"=c(NA,NA),
                    "p0.95"=c(NA,NA),
                    "p0.05"=c(NA,NA),
                    "p0.025"=c(NA,NA))
  return(c.d)
}

corTest <- function(d,s){
  
  # Create two datasets, one to add permutation columns to, one to pull from. because that's how I coded it?
  a <- d 
  perms <- d 
  nloci <- length(a$Locus)
  maxB <- max(a$BF)
  minB <- min(a$BF)
  maxG <- max(a$GL)
  minG <- min(a$GL)
  #print(length(a$Locus))
  
  # Safety check for if a subset of data has too few observations
  if(length(a$Locus) < 5){
    c.d <- emptyFrame(s)
    return(c.d)
  }
  
  # Define which column numbers to test for correlation
  
  b <- match('BF',names(a))
  g <- match('GL',names(a))
  
  # Permute BF values in 1,000 seperate columns, could permute BF or GL doesnt matter. 
  
  for (n in 1:1000) {
    newcolname1 <- paste('BF',n,sep='') 
    perms[, newcolname1] <- sample(a[[b]], replace = FALSE)
    #newcolname2 <- paste('GL',n,sep='') 
    #null[, newcolname2] <- sample(a[[g]], replace = FALSE)
  }
  
  # Create new dataset to add null rank coefficients to 
  
  assign('null', setNames(data.frame(1:1000), "replicate"))
  
  # Add bf gl 
  
  for (n in 1:1000) {
    
    # Grab column names, calculate value 
    colname1 <- paste('BF',n,sep='') 
    cc1.s <- cor(perms[,'GL'],perms[,colname1],method="spearman") 
    null[n,'null.spearman'] <- cc1.s
    cc1.p <- cor(perms[,'GL'],perms[,colname1],method="pearson") 
    null[n,'null.pearson'] <- cc1.p
  }
  
  
  # Dataframe for emp scores
  emp <- data.frame(TestType=character(),
                    corr.eff=integer(),
                    p=integer(),
                    p0.975=integer(),
                    p0.95=integer(),
                    p0.05=integer(),
                    p0.025=integer(),stringsAsFactors=FALSE)
  
  # Designate which row is for which test 
  
  emp[1,'TestType'] <- 'emp.spearman'
  emp[2,'TestType'] <- 'emp.pearson'
  
  
  # Calculate correlation stats from emperical data,  add to df
  
  e.s <- cor(a[[b]],a[[g]],method="spearman") 
  e.p <- cor(a[[b]],a[[g]],method="pearson") 
  
  emp[1,'corr.eff'] <- e.s
  emp[2,'corr.eff'] <- e.p
  
  
  # Calculate 2 sided p values
  
  p.s <- sum(abs(null$null.spearman) >= abs(e.s)) / length(null$null.spearman)
  p.p <- sum(abs(null$null.pearson) >= abs(e.s)) / length(null$null.pearson)
  
  emp[1,'p'] <- p.s
  emp[2,'p'] <- p.p
  
  
  # Add rows for quantile cut offs 
  
  emp[1,'p0.025'] <- quantile(null$null.spearman, c(0.025, 0.975),na.rm = TRUE)[[1]]
  emp[1,'p0.975'] <- quantile(null$null.spearman, c(0.025, 0.975),na.rm = TRUE)[[2]]
  emp[1,'p0.05'] <- quantile(null$null.spearman, c(0.05, 0.95),na.rm = TRUE)[[1]]
  emp[1,'p0.95'] <- quantile(null$null.spearman, c(0.05, 0.95),na.rm = TRUE)[[2]]
  emp[2,'p0.025'] <- quantile(null$null.pearson, c(0.025, 0.975),na.rm = TRUE)[[1]]
  emp[2,'p0.975'] <- quantile(null$null.pearson, c(0.025, 0.975),na.rm = TRUE)[[2]]
  emp[2,'p0.05'] <- quantile(null$null.pearson, c(0.05, 0.95),na.rm = TRUE)[[1]]
  emp[2,'p0.95'] <- quantile(null$null.pearson, c(0.05, 0.95),na.rm = TRUE)[[2]]
  
  # Add significance yes/no and other data 
  
  corr.data <- emp %>%
    mutate_if(is.numeric,round,2) %>% 
    add_column(Hypothesis = "ALL", .after = "TestType") %>% 
    add_column(Sample = s, .after = "TestType") %>% 
    add_column(Loci = nloci, .after = "p") %>%  
    add_column(minG = minG, .after = "Loci")%>% 
    add_column(maxG = maxG, .after = "Loci")%>% 
    add_column(minB = minB, .after = "Loci")%>% 
    add_column(maxB = maxB, .after = "Loci")
  
  return(corr.data)
}


setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/Tables")

bigCor <- corTest(all,"all") 


bigCor.0x <- corTest(all.0x,"sig0x")
bigCor.1x <- corTest(all.1x,"sig1x")
bigCor.2x <- corTest(all.2x,"sig2x")
bigCor.50m <- corTest(all.50m,"quant50m")
bigCor.60m <- corTest(all.60m,"quant60m")
bigCor.80m <- corTest(all.80m,"quant80m")
bigCor.50l <- corTest(all.50l,"50l")
bigCor.50h <- corTest(all.50h,"50h")




results <- bind_rows(bigCor,bigCor.0x, bigCor.1x, bigCor.2x, 
                     bigCor.50m, bigCor.60m, bigCor.80m, 
                     bigCor.50l, bigCor.50h) %>% 
  filter(TestType == "emp.spearman")


library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/Tables/AllSchmooshNoHypoth.xlsx",sep='')

write.xlsx2(as.data.frame(results), file=fname, row.names=FALSE)

