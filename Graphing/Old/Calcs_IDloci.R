rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library("viridis")
# set dataset 

#dataset <- "SinghalOG"
dataset <- "Streicher"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Colors 
# https://www.w3schools.com/colors/colors_picker.asp
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS"))
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))

# Get total number of loci 
loci <- length(mLGL$Locus)

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10 - not including average calcs
names(mLBF[7:10])
mLBF[7:10] <- mLBF[7:10]/2
# for singhal OG 
#names(mLBF[4])
#mLBF[4] <- mLBF[4]/2

# Create dataframe with all metadata, locus names, one col for hypothesis, one call for bf, and one gl
keep <- c(1,seq(14,22))

B <- bind_rows(reformatdf(mLBF,"BF","AIvSA",1,keep),
               reformatdf(mLBF,"BF","AIvSI",1,keep),
               reformatdf(mLBF,"BF","SAvSI",1,keep),
               reformatdf(mLBF,"BF","TvS",1,keep))

G <- bind_rows(reformatdf(mLGL,"GL","AIvSA",1,keep),
               reformatdf(mLGL,"GL","AIvSI",1,keep),
               reformatdf(mLGL,"GL","SAvSI",1,keep),
               reformatdf(mLGL,"GL","TvS",1,keep))
# for SinghalOG 
#keep <- c(1,seq(5,13)) 
#B <- reformatdf(mLBF,"BF","TvS",1,keep)
#G <- reformatdf(mLGL,"GL","TvS",1,keep)

metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")
all <- left_join(B, G, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))
some <- all %>% select(-metaData)


a <- some %>% mutate(zero = case_when(between(GL,-0.001,0.001) ~ 'z'))

a <- some %>% filter(between(GL,-0.001,0.001))

