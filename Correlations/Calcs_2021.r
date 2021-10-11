rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


# set dataset 

#dataset <- "Streicher"
#dataset <- "Singhal"
#dataset <- "Reeder"
dataset <- "Burbrink"
#dataset <- "Singhal"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data

## Read in data tables, drop any rows that have all NAs
mLBF <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/",dataset,"MargLikeSummaryNA64mil.txt",sep=""),
                   header=TRUE) %>% drop_na()
mLGL <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/",dataset,"MaxLikeSummaryNA.txt",sep=""),
                     header=TRUE) %>% drop_na()


# Check number of NA 
sum(is.na(mLGL$Sclero))

## Number of overlapping
length(intersect(mLBF$Locus,mLGL$Locus))

# Remove any other from the larger dataset - look in environ in R window to see who is bigger 
mLBF <- mLBF[mLBF$Locus %in% intersect(mLBF$Locus,mLGL$Locus),]
mLGL <- mLGL[mLGL$Locus %in% intersect(mLBF$Locus,mLGL$Locus),]


# Bayes factors 

# K = lnBF(m1,m2) = ln(m1) - ln(m2)
# 2ln(BF) = 2ln(e^K) = 2*K = 2(ln(m1) - ln(m2))

mL <- mLBF

## Calculate BF and put in new columns
mL$TvS <- 2*(mL$ToxPoly-mL$Sclero)
mL$AIvSA <- 2*(mL$ToxAngIg-mL$ToxSnAng)
mL$AIvSI <- 2*(mL$ToxAngIg-mL$ToxSnIg)
mL$SAvSI <- 2*(mL$ToxSnAng-mL$ToxSnIg)


# check NA 
sum(is.na(mL))

# Replace NA with zeros 

mL[is.na(mL)] <- 0
sum(is.na(mL))

## Switch back into mLBF and clear mL variable to use with dGLS

mLBF <- mL

rm(mL,all,n1,d1,x1,n2,d2,x2,n3,d3,x3)



# Different in gene log likelihood  

# dGLS = LH1 - LH2

mL <- mLGL

## Calculate dGLS and put in new columns
mL$TvS <- (mL$ToxPoly-mL$Sclero)
mL$AIvSA <- (mL$ToxAngIg-mL$ToxSnAng)
mL$AIvSI <- (mL$ToxAngIg-mL$ToxSnIg)
mL$SAvSI <- (mL$ToxSnAng-mL$ToxSnIg)
## Averages 
mL$AIvAvg <- (mL$ToxAngIg-((mL$ToxSnIg+mL$ToxSnAng)/2))
mL$SAvAvg <- (mL$ToxSnAng-((mL$ToxAngIg+mL$ToxSnIg)/2))
mL$SIvAvg <- (mL$ToxSnIg-((mL$ToxAngIg+mL$ToxSnAng)/2))

# check NA 
sum(is.na(mL))

# Replace NA with zeros 

mL[is.na(mL)] <- 0
sum(is.na(mL))


## Switch back into mLBF and clear mL variable to use with dGLS

mLGL <- mL

rm(mL)

# If you want to save..

save(mLGL, mLBF, file=paste("Calcs_",dataset,".RData",sep=""))
#save(mLBF, file=paste("Calcs_",dataset,"_BFonly.RData",sep=""))



