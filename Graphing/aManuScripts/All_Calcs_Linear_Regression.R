rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/Tables")


#y-gl
#x-bf
lr.df<- data.frame(hypothesis=character(),
                   subset=character(),
                   a=numeric(), 
                   b=numeric(), 
                   r2=numeric(),
                   nloci=numeric(),
                   source=character())
dsets <- list(all=all,all.2x=all.2x,all.1x=all.1x,all.0x=all.0x,
              all.80m=all.80m,all.60m=all.60m,all.80mb=all.80mb,all.60mb=all.60mb,all.50m=all.50m,
              all.50l=all.50l,all.50h=all.50h) 
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

fname <- "/Users/ChatNoir/Projects/Squam/Graphs/Tables/All_LinearRegression.xlsx"

write.xlsx2(as.data.frame(lr.df), file=fname, sheetName="LinearRegression", row.names=FALSE)

