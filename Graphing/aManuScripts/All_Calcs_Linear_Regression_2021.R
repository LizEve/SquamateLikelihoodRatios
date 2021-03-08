rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh2021.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/Tables")


# Need to redo subsets without splitting by dataset ------------------------------------------------------------------------------------------------------------------------------------------------------


# all bf within and all gl within. agnostic to other value. gl in OR bf in 
all.60m <- all %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) | 
           between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80m <- all %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) | 
           between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# both within - gl AND bf within  
all.60mb <- all %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) , 
         between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80mb <- all %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) , 
         between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

#y-gl
#x-bf
lr.df<- data.frame(hypothesis=character(),
                   subset=character(),
                   a=numeric(), 
                   b=numeric(), 
                   r2=numeric(),
                   nloci=numeric(),
                   source=character())
dsets <- list(all=all,
              all.80m=all.80m,all.60m=all.60m,all.80mb=all.80mb,all.60mb=all.60mb)

#dsets <- list(all=all,all.sx4=all.sx4,all.sx2=all.sx2) 

hypos <- c("TvS","AIvSA","AIvSI","SAvSI")
#hypos <- c("AIvSA","AIvSI","SAvSI") # reeder

for (x in seq(1,length(dsets))){
  dd <- names(dsets[x])
  d <- as.data.frame(dsets[x])
  v <- paste(dd,"variable",sep=".")
  #print(dd)
  #print(v)
  for (h in hypos){
    df <- d %>% filter(d[[v]] == h)
    b <- paste(dd,"BF",sep=".")
    g <- paste(dd,"GL",sep=".")
    m = lm(df[[b]] ~ df[[g]])
    print(dd)
    print(h)
    print(m)
    print(summary(m)$r.squared)
    l <- list(hypothesis=h, 
              subset=dd, 
              a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2), 
              r2 = format(summary(m)$r.squared, digits = 3),
              nloci = length(df[[g]]))
    #print(data.frame(l))
    lr.df <- rbind(lr.df,data.frame(l))
  }
}


library(xlsx)

fname <- "/Users/ChatNoir/Projects/Squam/Graphs/Tables/Combined_LinearRegression.xlsx"

write.xlsx2(as.data.frame(lr.df), file=fname, sheetName="LinearRegression", row.names=FALSE)

