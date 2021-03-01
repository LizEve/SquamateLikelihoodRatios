
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

dataset <- "Reeder"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")

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
all.short <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-20,20)) %>% 
  filter(between(GL,-8,8)) # "significance levels (5,0.5) x 4
all.shorter <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-10,10)) %>% 
  filter(between(GL,-8,8)) # "significance levels (5,0.5) x 2
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
    mutate(sig.level = case_when(corr.eff < p0.025 ~ "-0.05",
                                 corr.eff > p0.975 ~ "0.05",
                                 corr.eff >= p0.025 & corr.eff < p0.05 ~ "-0.10",
                                 corr.eff >= p0.95 & corr.eff < p0.975 ~ "0.10",
                                 corr.eff >= p0.05 & corr.eff <= p0.95 ~ "None")) %>%
    mutate_if(is.numeric,round,2) %>% 
    add_column(Hypothesis = h, .after = "TestType") %>% 
    add_column(Sample = s, .after = "TestType") %>% 
    add_column(Loci = nloci)
  return(corr.data)
}


tvs <- corTest(all.exp,"TvS","all") 
aivsa <- corTest(all.exp,"AIvSA","all")
aivsi <- corTest(all.exp,"AIvSI","all")
savsi <- corTest(all.exp,"SAvSI","all")

tvs.ser <- corTest(all.shorter,"TvS","sigx2")
aivsa.ser <- corTest(all.shorter,"AIvSA","sigx2")
aivsi.ser <- corTest(all.shorter,"AIvSI","sigx2")
savsi.ser <- corTest(all.shorter,"SAvSI","sigx2")

tvs.s <- corTest(all.short,"TvS","sigx4")
aivsa.s <- corTest(all.short,"AIvSA","sigx4")
aivsi.s <- corTest(all.short,"AIvSI","sigx4")
savsi.s <- corTest(all.short,"SAvSI","sigx4")

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
                     tvs.ser,aivsa.ser,aivsi.ser,savsi.ser,
                     tvs.s,aivsa.s,aivsi.s,savsi.s) %>% filter(TestType == "emp.spearman")

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

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_correlation2.xlsx",sep='')

write.xlsx2(as.data.frame(results), file=fname, sheetName="ALL", row.names=FALSE)
