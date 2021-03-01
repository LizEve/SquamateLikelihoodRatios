rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


BFofAvg <- function(df,numor,denoms){
  
  # Create a new dataset for new BF values & add a new empty column to populate later 
  new <- subset(df, select = c('Locus'))
  new$newestAvgBF <- vector(mode="numeric",length=length(new$Locus))
  
  # for locus in dataset 
  for (i in 1:nrow(new)){
    
    # Vector of denoms for this i 
    idenoms <- denoms[i,c(2,length(names(denoms)))]
    
    # Maximum of alternate margLH
    a <- max(idenoms)
    
    # Calculate average across alternate constraints - sum does a sigma sum across all values in the list passed to it
    avg <- a+log(sum(exp(idenoms-a)))-log(length(idenoms))
    
    # Calculate bayes factor - 2ln(BF)
    new[i,2] <- 2*(numor[i,2] - avg)
    }
  
  return(new)
}



# set dataset 

#dataset <- "Streicher"
#dataset <- "Singhal"
#dataset <- "Reeder"
dataset <- "Burbrink"
#dataset <- "SinghalOG"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data

## Read in data tables, drop any rows that have all NAs
mLBF <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/",dataset,"MargLikeSummaryNA.txt",sep=""),
                   header=TRUE) %>% drop_na()
mLGL <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/",dataset,"MaxLikeSummaryNA.txt",sep=""),
                     header=TRUE) %>% drop_na()
#mLGL <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/",dataset,"MaxLikeSummary.txt",sep=""),
#                   header=TRUE) %>% drop_na()


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
mL$SAvSI <- 2*(mL$ToxSnAng-mL$ToxSnIg)

# Select numorator and denominators to average over 
# Get BF values 
# Add new values to main dataframe 
n1 <- mL[,(names(mL) %in% c("ToxAngIg","Locus"))]
d1 <- mL[,(names(mL) %in% c("ToxSnAng", "ToxSnIg","Locus"))]
x1 <- BFofAvg(mL,n1,d1)
mL<-left_join(mL,x1,by='Locus') %>% rename(AIvAvg=newestAvgBF)

# Select numorator and denominators to average over 
# Get BF values 
# Add new values to main dataframe 
n2 <- mL[,(names(mL) %in% c("ToxSnAng","Locus"))]
d2 <- mL[,(names(mL) %in% c("ToxAngIg", "ToxSnIg","Locus"))]
x2 <- BFofAvg(mL,n2,d2)
mL<-left_join(mL,x2,by='Locus') %>% rename(SAvAvg=newestAvgBF)

# Select numorator and denominators to average over 
# Get BF values 
# Add new values to main dataframe 
n3 <- mL[,(names(mL) %in% c("ToxSnIg","Locus"))]
d3 <- mL[,(names(mL) %in% c("ToxSnAng", "ToxAngIg","Locus"))]
x3 <- BFofAvg(mL,n3,d3)
mL<-left_join(mL,x3,by='Locus') %>% rename(SIvAvg=newestAvgBF)


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
mL$SAvSI <- (mL$ToxSnAng-mL$ToxSnIg)
## Averages 
mL$AIvAvg <- (mL$ToxAngIg-((mL$ToxSnIg+mL$ToxSnAng)/2))
mL$SAvAvg <- (mL$ToxSnAng-((mL$ToxAngIg+mL$ToxSnIg)/2))
mL$SIvAvg <- (mL$ToxSnIg-((mL$ToxAngIg+mL$ToxSnAng)/2))

## Switch back into mLBF and clear mL variable to use with dGLS

mLGL <- mL

rm(mL)


# If you want to save..

save(mLGL, mLBF, file=paste("Calcs_",dataset,".RData",sep=""))
#save(mLBF, file=paste("Calcs_",dataset,"_BFonly.RData",sep=""))



