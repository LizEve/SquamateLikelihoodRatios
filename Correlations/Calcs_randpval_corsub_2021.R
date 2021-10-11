

rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(readxl)

setwd("/Users/ChatNoir/Projects/Squam/Graphs/")

# Read in data
dfmid <- read_excel("/Users/chatnoir/Projects/Squam/Graphs/Tables/Correlation_forR_2021.xlsx")
df <- dfmid %>%  mutate(Sample = sub('sig', '', Sample))
mdf <- data.frame(Sample=character(0),Hypothesis=character(),corr.eff=numeric(0),pRand=numeric(0),dataset=character(0))
#------------------------------------------------------------------------------------------------------------------------------------------------------

# Streicher ---------------------------------------------- Streicher Pvals --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsampleStreicher")
results <- results %>% filter(Sample != "all")

# Remove non subsampled values 
d <- "Streicher"
mid <- df %>% filter(Sample != "all") %>% filter(dataSet == d)

hyps <- c('TvS','AIvSA', 'AIvSI','SAvSI')
subs <- c('80m','60m','80mb','60mb')
l1 <- 300
l2 <- 2000

for (h in hyps) {
  for (s in subs){
    if (h == 'TvS'){
      l <- l1
    }else{
      l <- l2}
    n <- results %>% filter(Hypothesis == h) %>% filter(Loci == l)
    m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
    e <- m$corr.eff
    p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
    a <- data.frame(Sample=s,Hypothesis=h,corr.eff=e,pRand=p.s,dataset=d)
    mdf <- rbind(mdf,a)
  }
}


# all p = 1, ie all sim cor eff are higher than emp vals. all of them. 



#------------------------------------------------------------------------------------------------------------------------------------------------------

# Singhal ---------------------------------------------- Singhal --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsampleSinghal")
results <- results %>% filter(Sample != "all")

# Remove non subsampled values 
d <- "Singhal"
mid <- df %>% filter(Sample != "all") %>% filter(dataSet == d)

hyps <- c('TvS','AIvSA', 'AIvSI','SAvSI')
subs <- c('80m','60m','80mb','60mb')
l1 <- 1000
l2 <- 3000

for (h in hyps) {
  for (s in subs){
    if (h == 'TvS'){
      l <- l1
    }else{
      l <- l2}
    n <- results %>% filter(Hypothesis == h) %>% filter(Loci == l)
    m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
    e <- m$corr.eff
    p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
    a <- data.frame(Sample=s,Hypothesis=h,corr.eff=e,pRand=p.s,dataset=d)
    mdf <- rbind(mdf,a)
  }
}

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Reeder ---------------------------------------------- Reeder --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsampleReeder")
results <- results %>% filter(Sample != "all")

# Remove non subsampled values 
d <- "Reeder"
mid <- df %>% filter(Sample != "all") %>% filter(dataSet == d)

hyps <- c('TvS','AIvSA', 'AIvSI','SAvSI')
subs <- c('80m','60m','80mb','60mb')
l1 <- 30
l2 <- 25

for (h in hyps) {
  for (s in subs){
    if (h == 'TvS'){
      l <- l1
    }else{
      l <- l2}
    n <- results %>% filter(Hypothesis == h) %>% filter(Loci == l)
    m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
    e <- m$corr.eff
    p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
    a <- data.frame(Sample=s,Hypothesis=h,corr.eff=e,pRand=p.s,dataset=d)
    mdf <- rbind(mdf,a)
  }
}


#------------------------------------------------------------------------------------------------------------------------------------------------------

# Burbrink ---------------------------------------------- Burbrink --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsampleBurbrink")
results <- results %>% filter(Sample != "all")

# Remove non subsampled values 
d <- "Burbrink"
mid <- df %>% filter(Sample != "all") %>% filter(dataSet == d)

hyps <- c('TvS','AIvSA', 'AIvSI','SAvSI')
subs <- c('80m','60m','80mb','60mb')
l1 <- 100
l2 <- 200

for (h in hyps) {
  for (s in subs){
    if (h == 'TvS'){
      l <- l1
    }else{
      l <- l2}
    n <- results %>% filter(Hypothesis == h) %>% filter(Loci == l)
    m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
    e <- m$corr.eff
    p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
    a <- data.frame(Sample=s,Hypothesis=h,corr.eff=e,pRand=p.s,dataset=d)
    mdf <- rbind(mdf,a)
  }
}

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Save ---------------------------------------------- Save --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------


library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/Tables/Corr_rand_pvals_mmb_2021.xlsx",sep='')

write.xlsx2(as.data.frame(mdf), file=fname, row.names=FALSE)


