
rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)

# set dataset 

## In manuscript and edited xlsx sheets i refer to 1x as 2x 

dataset <- "Singhal"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS")) %>%
  mutate_if(is.numeric,round,2)
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF")) %>%
  mutate_if(is.numeric,round,2)

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10
mLBF[7:10] <- mLBF[7:10]/2

# Transform datasets - stack two comparisons, and add column for bf and gl 
# Each dataset one hypothesis 

B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))


G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

# Smoosh together again
all <- bind_rows(B,G) 



# Originally ran with threshold of dGLS being 0.5, 2x = 1, 4x=2
# trying again with 2, 2x=4, 4x=8

all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
all.exp <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% select(-c("supportType.b","supportType.g"))


all.2x <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-15,15) | 
           between(GL,-6,6)) # "significance levels (5,2) 5 + 5x 2
all.1x <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-10,10) | 
           between(GL,-4,4)) # "significance levels (5,2) 5 + 5x 1
all.0x <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-5,5) | 
           between(GL,-2,2)) # "significance levels (5,2) 5 + 5x 1


# get percentiles 

all.50l <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0),quantile(all.exp$GL, 0.5)) |
           between(BF,quantile(all.exp$BF,0),quantile(all.exp$BF, 0.5)))

all.50h <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0.5),quantile(all.exp$GL, 1)) |
           between(BF,quantile(all.exp$BF,0.5),quantile(all.exp$BF, 1)))

all.50m <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0.25),quantile(all.exp$GL, 0.75)) |
           between(BF,quantile(all.exp$BF,0.25),quantile(all.exp$BF, 0.75)))

all.60m <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0.2),quantile(all.exp$GL, 0.8)) | 
           between(BF,quantile(all.exp$BF,0.2),quantile(all.exp$BF, 0.8)))

all.80m <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0.1),quantile(all.exp$GL, 0.9)) | 
           between(BF,quantile(all.exp$BF,0.1),quantile(all.exp$BF, 0.9)))
  

#  filter(between(GL,quantile(all.exp$GL, 0.25),quantile(all.exp$GL, 0.75))) # Middle 50% of distribution



# Filter based on only one threshold
#all.short.b <- all.exp %>% group_by(variable) %>% 
#  filter(between(BF,-20,20)) # "significance levels (5,0.5) x 4
#all.shorter.b <- all.exp %>% group_by(variable) %>% 
#  filter(between(BF,-10,10))  # "significance levels (5,0.5) x 2
#all.short.g <- all.exp %>% group_by(variable) %>% 
#  filter(between(GL,-2,2)) # "significance levels (5,0.5) x 4
#all.shorter.g <- all.exp %>% group_by(variable) %>% 
#  filter(between(GL,-4,4)) # "significance levels (5,0.5) x 2
#all.mid.b <- all.exp %>% group_by(variable) %>% 
#  filter(between(BF,quantile(all.exp$BF, 0.25),quantile(all.exp$BF, 0.75)))  # Middle 50% of distribution
#all.mid.g <- all.exp %>% group_by(variable) %>% 
#  filter(between(GL,quantile(all.exp$GL, 0.25),quantile(all.exp$GL, 0.75))) # Middle 50% of distribution

# Testing total number of loci - 1/3 of total loci
# filter(variable != 'TvS') %>% 
#all.exp <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% select(-c("supportType.b","supportType.g"))
#all.short <- all.exp %>% group_by(variable) %>% 
#  filter(between(BF,-20,20)) %>% 
#  filter(between(GL,-2,2)) # "significance levels (5,0.5) x 4
#all.shorter <- all.exp %>% group_by(variable) %>% 
#  filter(between(BF,-10,10)) %>% 
#  filter(between(GL,-1,1)) # "significance levels (5,0.5) x 2
#nl <- 30
#all.exp <- all.exp %>% group_by(variable) %>%  sample_n(nl)
#all.short <- all.short %>%  group_by(variable) %>%sample_n(nl)
#all.shorter <- all.shorter %>% group_by(variable) %>%sample_n(nl)
#
##proof of concept that Im getting the mid two quantiles
#h <- "AIvSA"
#datf <- all.exp[all.exp$variable==h,]
#d.b <- datf$BF
#d.g <- datf$GL
#b <- IQR.outliers(datf,d.b)
#g <- IQR.outliers(datf,d.g)
#quantile(d.b)[1]
#quantile(d.b)[5]
#quantile(d.b)[2]
#quantile(d.b)[4]
#quantile(datf$BF, 0.25)
#quantile(datf$BF, 0.75)

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

tvs <- corTest(all.exp,"TvS","all") 
aivsa <- corTest(all.exp,"AIvSA","all")
aivsi <- corTest(all.exp,"AIvSI","all")
savsi <- corTest(all.exp,"SAvSI","all")

tvs.1x <- corTest(all.1x,"TvS","sig1x")
aivsa.1x <- corTest(all.1x,"AIvSA","sig1x")
aivsi.1x <- corTest(all.1x,"AIvSI","sig1x")
savsi.1x <- corTest(all.1x,"SAvSI","sig1x")

tvs.2x <- corTest(all.2x,"TvS","sig2x")
aivsa.2x <- corTest(all.2x,"AIvSA","sig2x")
aivsi.2x <- corTest(all.2x,"AIvSI","sig2x")
savsi.2x <- corTest(all.2x,"SAvSI","sig2x")


tvs.50m <- corTest(all.50m,"TvS","sig50m")
aivsa.50m <- corTest(all.50m,"AIvSA","sig50m")
aivsi.50m <- corTest(all.50m,"AIvSI","sig50m")
savsi.50m <- corTest(all.50m,"SAvSI","sig50m")


tvs.60m <- corTest(all.60m,"TvS","sig60m")
aivsa.60m <- corTest(all.60m,"AIvSA","sig60m")
aivsi.60m <- corTest(all.60m,"AIvSI","sig60m")
savsi.60m <- corTest(all.60m,"SAvSI","sig60m")


tvs.80m <- corTest(all.80m,"TvS","sig80m")
aivsa.80m <- corTest(all.80m,"AIvSA","sig80m")
aivsi.80m <- corTest(all.80m,"AIvSI","sig80m")
savsi.80m <- corTest(all.80m,"SAvSI","sig80m")


tvs.0x <- corTest(all.0x,"TvS","0x")
aivsa.0x <- corTest(all.0x,"AIvSA","0x")
aivsi.0x <- corTest(all.0x,"AIvSI","0x")
savsi.0x <- corTest(all.0x,"SAvSI","0x")


tvs.50l <- corTest(all.50l,"TvS","50l")
aivsa.50l <- corTest(all.50l,"AIvSA","50l")
aivsi.50l <- corTest(all.50l,"AIvSI","50l")
savsi.50l <- corTest(all.50l,"SAvSI","50l")

tvs.50h <- corTest(all.50h,"TvS","50h")
aivsa.50h <- corTest(all.50h,"AIvSA","50h")
aivsi.50h <- corTest(all.50h,"AIvSI","50h")
savsi.50h <- corTest(all.50h,"SAvSI","50h")



## Filter by only one value threshold
#tvs.ser.b <- corTest(all.shorter.b,"TvS","sigx2.b")
#aivsa.ser.b <- corTest(all.shorter.b,"AIvSA","sigx2.b")
#aivsi.ser.b <- corTest(all.shorter.b,"AIvSI","sigx2.b")
#savsi.ser.b <- corTest(all.shorter.b,"SAvSI","sigx2.b")
#
#tvs.s.b <- corTest(all.short.b,"TvS","sigx4.b")
#aivsa.s.b <- corTest(all.short.b,"AIvSA","sigx4.b")
#aivsi.s.b <- corTest(all.short.b,"AIvSI","sigx4.b")
#savsi.s.b <- corTest(all.short.b,"SAvSI","sigx4.b")
#
#tvs.ser.g <- corTest(all.shorter.g,"TvS","sigx2.g")
#aivsa.ser.g <- corTest(all.shorter.g,"AIvSA","sigx2.g")
#aivsi.ser.g <- corTest(all.shorter.g,"AIvSI","sigx2.g")
#savsi.ser.g <- corTest(all.shorter.g,"SAvSI","sigx2.g")
#
#tvs.s.g <- corTest(all.short.g,"TvS","sigx4.g")
#aivsa.s.g <- corTest(all.short.g,"AIvSA","sigx4.g")
#aivsi.s.g <- corTest(all.short.g,"AIvSI","sigx4.g")
#savsi.s.g <- corTest(all.short.g,"SAvSI","sigx4.g")
#
#tvs.m.b <- corTest(all.mid.b,"TvS","mid50.b")
#aivsa.m.b <- corTest(all.mid.b,"AIvSA","mid50.b")
#aivsi.m.b <- corTest(all.mid.b,"AIvSI","mid50.b")
#savsi.m.b <- corTest(all.mid.b,"SAvSI","mid50.b")
#
#tvs.m.g <- corTest(all.mid.g,"TvS","mid50.g")
#aivsa.m.g <- corTest(all.mid.g,"AIvSA","mid50.g")
#aivsi.m.g <- corTest(all.mid.g,"AIvSI","mid50.g")
#savsi.m.g <- corTest(all.mid.g,"SAvSI","mid50.g")


results <- bind_rows(tvs,aivsa,aivsi,savsi,
                     tvs.80m,aivsa.80m,aivsi.80m,savsi.80m,
                     tvs.60m,aivsa.60m,aivsi.60m,savsi.60m,
                     tvs.50m,aivsa.50m,aivsi.50m,savsi.50m,
                     tvs.50l,aivsa.50l,aivsi.50l,savsi.50l,
                     tvs.50h,aivsa.50h,aivsi.50h,savsi.50h,
                     tvs.2x,aivsa.2x,aivsi.2x,savsi.2x,
                     tvs.1x,aivsa.1x,aivsi.1x,savsi.1x,
                     tvs.0x,aivsa.0x,aivsi.0x,savsi.0x) %>% filter(TestType == "emp.spearman")


##  Filter by only one value threshold
#results <- bind_rows(tvs,aivsa,aivsi,savsi,
#                     tvs.ser,aivsa.ser,aivsi.ser,savsi.ser,
#                     tvs.s,aivsa.s,aivsi.s,savsi.s,
#                     tvs.ser.b,aivsa.ser.b,aivsi.ser.b,savsi.ser.b,
#                     tvs.s.b,aivsa.s.b,aivsi.s.b,savsi.s.b,
#                     tvs.ser.g,aivsa.ser.g,aivsi.ser.g,savsi.ser.g,
#                     tvs.s.g,aivsa.s.g,aivsi.s.g,savsi.s.g,
#                     tvs.m.g,aivsa.m.g,aivsi.m.g,savsi.m.g,
#                     tvs.m.b,aivsa.m.b,aivsi.m.b,savsi.m.b) %>% filter(TestType == "emp.spearman")
#
#
#results <- bind_rows(tvs,aivsa,aivsi,savsi,
#                     tvs.ser,aivsa.ser,aivsi.ser,savsi.ser,
#                     tvs.s,aivsa.s,aivsi.s,savsi.s) %>% filter(TestType == "emp.spearman")

#results30a <- results
#results30b <- results
#results30c <- results

#results <- bind_rows(results30a, results30b, results30c)

#library(xlsx)

#fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_correlationSubsampled3.xlsx",sep='')

#write.xlsx2(as.data.frame(results), file=fname, sheetName="Subsample", row.names=FALSE)


# output excel sheet

library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_correlation_redo.xlsx",sep='')

write.xlsx2(as.data.frame(results), file=fname, sheetName="ALL", row.names=FALSE)

