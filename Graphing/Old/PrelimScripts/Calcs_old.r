rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# set dataset 

#dataset <- "Streicher"
#dataset <- "Singhal"
dataset <- "Reeder"
#dataset <- "Burbrink"


setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data

## Read in data tables, drop any rows that have all NAs
mLBF <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/",dataset,"MargLikeSummaryNA.txt",sep=""),
                   header=TRUE) %>% drop_na()
mLGL <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/",dataset,"MaxLikeSummaryNA.txt",sep=""),
                     header=TRUE) %>% drop_na()

## Number of overlapping
length(intersect(mLBF$Locus,mLGL$Locus))

# Remove any other from the larger dataset
mLBF <- mLBF[mLBF$Locus %in% intersect(mLBF$Locus,mLGL$Locus),]


# Bayes factors 

# K = lnBF(m1,m2) = ln(m1) - ln(m2)
# 2ln(BF) = 2ln(e^K) = 2*K = 2(ln(m1) - ln(m2))

mL <- mLBF

## Calculate BF and put in new columns
mL$TvS <- 2*(mL$ToxPoly-mL$Sclero)
mL$AIvSA <- 2*(mL$ToxAngIg-mL$ToxSnAng)
mL$AIvSI <- 2*(mL$ToxAngIg-mL$ToxSnIg)
mL$SAvAI <- 2*(mL$ToxSnAng-mL$ToxAngIg)
mL$SAvSI <- 2*(mL$ToxSnAng-mL$ToxSnIg)
mL$SIvAI <- 2*(mL$ToxSnIg-mL$ToxAngIg)
mL$SIvSA <- 2*(mL$ToxSnIg-mL$ToxSnAng)

## Make column of support interpretation for each hypothesis 

mL <- mutate(mL, TvS_support = case_when(TvS > 10 ~ "Very Strong for Tox",
                                         TvS <= 10 &  TvS >= 6 ~ "Strong for Tox",
                                         TvS > -6 & TvS < 6 ~ "Ambiguous",
                                         TvS >= -10 & TvS <= 6 ~ "Strong for Sclero",
                                         TvS < -10 ~ "Very Strong for Sclero"))

TvS_BFsupport <- mL %>% group_by(TvS_support) %>% dplyr::summarise(counts=n())
TvS_BFsupport$TvS_support <- factor(TvS_BFsupport$TvS_support, levels=c("Very Strong for Tox","Strong for Tox","Ambiguous","Strong for Sclero","Very Strong for Sclero"))

## Switch back into mLBF and clear mL variable to use with dGLS

mLBF <- mL

rm(mL)

# Different in gene log likelihood  

# dGLS = LH1 - LH2

mL <- mLGL

## Calculate dGLS and put in new columns
mL$TvS <- (mL$ToxPoly-mL$Sclero)
mL$AIvSA <- (mL$ToxAngIg-mL$ToxSnAng)
mL$AIvSI <- (mL$ToxAngIg-mL$ToxSnIg)
mL$SAvAI <- (mL$ToxSnAng-mL$ToxAngIg)
mL$SAvSI <- (mL$ToxSnAng-mL$ToxSnIg)
mL$SIvAI <- (mL$ToxSnIg-mL$ToxAngIg)
mL$SIvSA <- (mL$ToxSnIg-mL$ToxSnAng)

## Make column of support interpretation for each hypothesis 

mL <- mutate(mL, TvS_support = case_when(TvS > 0.5 ~ "Very Strong for Tox",
                                         TvS <= 0.5 &  TvS >= 0.1 ~ "Strong for Tox",
                                         TvS > -0.1 & TvS < 0.1 ~ "Ambiguous",
                                         TvS >= -0.5 & TvS <= 0.1 ~ "Strong for Sclero",
                                         TvS < -0.5 ~ "Very Strong for Sclero"))

TvS_GLsupport <- mL %>% group_by(TvS_support) %>% dplyr::summarise(counts=n())
TvS_GLsupport$TvS_support <- factor(TvS_GLsupport$TvS_support, levels=c("Very Strong for Tox","Strong for Tox","Ambiguous","Strong for Sclero","Very Strong for Sclero"))

## Switch back into mLBF and clear mL variable to use with dGLS

mLGL <- mL

rm(mL)


# If you want to save..

#save(mLGL, mLBF, TvS_GLsupport, TvS_BFsupport, file=paste("Calcs_",dataset,".RData",sep=""))
#save(mLBF, TvS_BFsupport, file=paste("Calcs_",dataset,"_BFonly.RData",sep=""))
