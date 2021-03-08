rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(readxl)
library(stats)
library(tibble)
library(reshape)
library(xlsx)

setwd("/Users/ChatNoir/Projects/Squam/Graphs/")

# Read in data
dfmid <- read_excel("/Users/chatnoir/Projects/Squam/Graphs/Tables/CombinedSchmoosh_2021.xlsx")
df <- dfmid %>%  mutate(Sample = sub('quant', '', Sample))
mdf <- data.frame(Sample=character(0),Hypothesis=character(),corr.eff=numeric(0),pRand=numeric(0))




load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsample")


results <- results %>% filter(Sample != "all")

# Remove non subsampled values 

mid <- df %>% filter(Sample != "all") 

hyps <- c('TvS','AIvSA', 'AIvSI','SAvSI')
subs <- c('80m','60m','80mb','60mb')
l1 <- 2000
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
    a <- data.frame(Sample=s,Hypothesis=h,corr.eff=e,pRand=p.s)
    mdf <- rbind(mdf,a)
  }
}


#------------------------------------------------------------------------------------------------------------------------------------------------------

# Save ---------------------------------------------- Save --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------


library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/Tables/Combined_corr_rand_pvals_mmb_2021.xlsx",sep='')

write.xlsx2(as.data.frame(mdf), file=fname, row.names=FALSE)




