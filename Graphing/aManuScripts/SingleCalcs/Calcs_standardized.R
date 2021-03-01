rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(moments)
library(reshape)

# set dataset 

dataset <- "Burbrink"

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
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
a <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% select(-c("supportType.b","supportType.g"))


# Z score - did this with dplyr and group in the graphs script
TvS <- a[a$variable=="TvS",]
TvS$BF.s <- scale(TvS$BF, center=TRUE, scale=TRUE)
TvS$GL.s <- scale(TvS$GL, center=TRUE, scale=TRUE)

AIvSA <- a[a$variable=="AIvSA",]
AIvSA$BF.s <- scale(AIvSA$BF, center=TRUE, scale=TRUE)
AIvSA$GL.s <- scale(AIvSA$GL, center=TRUE, scale=TRUE)

AIvSI <- a[a$variable=="AIvSI",]
AIvSI$BF.s <- scale(AIvSI$BF, center=TRUE, scale=TRUE)
AIvSI$GL.s <- scale(AIvSI$GL, center=TRUE, scale=TRUE)

SAvSI <- a[a$variable=="SAvSI",]
SAvSI$BF.s <- scale(SAvSI$BF, center=TRUE, scale=TRUE)
SAvSI$GL.s <- scale(SAvSI$GL, center=TRUE, scale=TRUE)

# Smoosh together again
all <- bind_rows(TvS,AIvSA,AIvSI,SAvSI) 

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Mean Mode Variance -STANDARDIZED-----------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

bf <- all %>% 
  group_by(variable) %>% 
  mutate( z.outlier = case_when(between(BF.s, mean(BF.s)-sd(BF.s)*3, mean(BF.s)+sd(BF.s)*3) ~ FALSE,
                                BF.s < mean(BF.s)-sd(BF.s)*3 ~ TRUE, 
                                BF.s > mean(BF.s)+sd(BF.s)*3 ~ TRUE)) %>%
  dplyr::summarise (Mean=mean(BF.s), Median=median(BF.s), Std=sd(BF.s), 
             Mad=mad(BF.s), IQR=IQR(BF.s), Max=max(BF.s), Min=min(BF.s),
             Outliers=sum(z.outlier),Out.Prop=sum(z.outlier)/n())

gl <- all %>% group_by(variable) %>% 
  mutate( z.outlier = case_when(between(GL.s, mean(GL.s)-sd(GL.s)*3, mean(GL.s)+sd(GL.s)*3) ~ FALSE,
                                GL.s < mean(GL.s)-sd(GL.s)*3 ~ TRUE, 
                                GL.s > mean(GL.s)+sd(GL.s)*3 ~ TRUE)) %>%
  dplyr::summarise(Mean=mean(GL.s), Median=median(GL.s), Std=sd(GL.s), 
            Mad=mad(GL.s), IQR=IQR(GL.s), Max=max(GL.s), Min=min(GL.s),
            Out.Sum=sum(z.outlier), Out.Prop=sum(z.outlier)/n()) 


stats.all <- all %>% 
  group_by(variable) %>% 
  mutate( z.b.outlier = case_when(between(BF.s, mean(BF.s)-sd(BF.s)*3, mean(BF.s)+sd(BF.s)*3) ~ FALSE,
                                BF.s < mean(BF.s)-sd(BF.s)*3 ~ TRUE, 
                                BF.s > mean(BF.s)+sd(BF.s)*3 ~ TRUE)) %>%
  mutate( z.g.outlier = case_when(between(GL.s, mean(GL.s)-sd(GL.s)*3, mean(GL.s)+sd(GL.s)*3) ~ FALSE,
                                GL.s < mean(GL.s)-sd(GL.s)*3 ~ TRUE, 
                                GL.s > mean(GL.s)+sd(GL.s)*3 ~ TRUE)) %>%
  dplyr::summarise(Total.Loci=n(),Mean.B=mean(BF.s), Mean.G=mean(GL.s), 
             Median.B=median(BF.s), Median.G=median(GL.s), 
             Std.B=sd(BF.s), Std.G=sd(GL.s),
             Max.B=max(BF.s), Max.G=max(GL.s),
             Min.B=min(BF.s), Min.G=min(GL.s),
             Out.Sum.B=sum(z.b.outlier), Out.Sum.G=sum(z.g.outlier),
             Out.Prop.B=sum(z.b.outlier)/n(), Out.Prop.G=sum(z.g.outlier)/n(),
             Mad.B=mad(BF.s), Mad.G=mad(GL.s), 
             IQR.B=IQR(BF.s), IQR.G=IQR(GL.s)) 
#%>% mutate_if(is.numeric,round,20)
# IQR = Q3-Q1 how spread out the middle of our data is

#Outliers based on IQR--STANDARDIZED-----------------------------------------------------------------------------------------------------------------------------------------------------
# Had to alter from <= to < because IQR for GL is 0 in some cases
# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r
# https://www.r-bloggers.com/combined-outlier-detection-with-dplyr-and-ruler/
require("dplyr")
IQR.outliers <- function(df,d){
  lowerq = quantile(d)[2]
  upperq = quantile(d)[4]
  iqr = upperq - lowerq 
  extreme.threshold.upper = (iqr * 1.5) + upperq
  extreme.threshold.lower = lowerq - (iqr * 1.5)
  z <- df %>% mutate(outlier = case_when(d  < extreme.threshold.lower ~ TRUE,
                                    d > extreme.threshold.upper ~ TRUE,
                                    between(d,extreme.threshold.lower,extreme.threshold.upper) ~ FALSE)) %>% 
    dplyr::group_by(outlier) %>% dplyr::count()
  return(z)
}


outlier.counter <- function(h){
  datf <- all[all$variable==h,]
  d.b <- datf$BF.s
  d.g <- datf$GL.s
  b <- IQR.outliers(datf,d.b)
  g <- IQR.outliers(datf,d.g)
  return(c(b,g))
}

# Prints number of outliers according to IQR, bf then gl 
# AIvSA, AIvSI, SAvSI, TvS
stats.all$Out.Iqr.B <- c(outlier.counter('AIvSA')[2][[1]][2],
                         outlier.counter('AIvSI')[2][[1]][2],
                         outlier.counter('SAvSI')[2][[1]][2],
                         outlier.counter('TvS')[2][[1]][2])

stats.all$Out.Iqr.G <- c(outlier.counter('AIvSA')[4][[1]][2],
                         outlier.counter('AIvSI')[4][[1]][2],
                         outlier.counter('SAvSI')[4][[1]][2],
                         outlier.counter('TvS')[4][[1]][2])
# Get proportion of loci 
stats.all$Out.Iqr.Prop.B <- stats.all$Out.Iqr.B/stats.all$Total.Loci
stats.all$Out.Iqr.Prop.G <- stats.all$Out.Iqr.G/stats.all$Total.Loci

stats.all$variable =paste(stats.all$variable, ".std",sep='')

#---------------------------------------------------------------------------------------------------------------------------------------------------
# NON STANDARD------------------------------------------------------------------------------------------------------------------------------------------------------

# Mean Mode Variance Skew - NON STANDARD-------------------------------------------------------------------------------------------------------------------------------------------------
#Skew and kurtosis same for standard and original values
#https://brownmath.com/stat/shape.htm#SkewnessInterpret

# > 3, leptokurtic, sharp peak in graph
# = 3, mesokurtic, normal distribution
# < 3 platykurtic, flatter/wider curve

# positive - skew right, longer right tail
# negative - skew left, left tail longer 
# > +-1 - greatly skewed 
# 0.5 - 1 - moderately skewed
# < +-0.5 - approx symmetric 

bf.n <- all %>% 
  group_by(variable) %>% 
  mutate( z.outlier = case_when(between(BF, mean(BF)-sd(BF)*3, mean(BF)+sd(BF)*3) ~ FALSE,
                                BF < mean(BF)-sd(BF)*3 ~ TRUE, 
                                BF > mean(BF)+sd(BF)*3 ~ TRUE)) %>%
  dplyr::summarise (Mean=mean(BF), Median=median(BF), Std=sd(BF), 
             Mad=mad(BF), IQR=IQR(BF), Max=max(BF), Min=min(BF),
             Outliers=sum(z.outlier),Out.Prop=sum(z.outlier)/n(), 
             Skew = skewness(BF), Kurt = kurtosis(BF))

gl.n <- all %>% group_by(variable) %>% 
  mutate( z.outlier = case_when(between(GL, mean(GL)-sd(GL)*3, mean(GL)+sd(GL)*3) ~ FALSE,
                                GL < mean(GL)-sd(GL)*3 ~ TRUE, 
                                GL > mean(GL)+sd(GL)*3 ~ TRUE)) %>%
  dplyr::summarise(Mean=mean(GL), Median=median(GL), Std=sd(GL), 
            Mad=mad(GL), IQR=IQR(GL), Max=max(GL), Min=min(GL),
            Out.Sum=sum(z.outlier), Out.Prop=sum(z.outlier)/n(),
            Skew = skewness(GL), Kurt = kurtosis(GL)) 


stats.all.n <- all %>% 
  group_by(variable) %>% 
  mutate( z.b.outlier = case_when(between(BF, mean(BF)-sd(BF)*3, mean(BF)+sd(BF)*3) ~ FALSE,
                                  BF < mean(BF)-sd(BF)*3 ~ TRUE, 
                                  BF > mean(BF)+sd(BF)*3 ~ TRUE)) %>%
  mutate( z.g.outlier = case_when(between(GL, mean(GL)-sd(GL)*3, mean(GL)+sd(GL)*3) ~ FALSE,
                                  GL < mean(GL)-sd(GL)*3 ~ TRUE, 
                                  GL > mean(GL)+sd(GL)*3 ~ TRUE)) %>%
  dplyr::summarise(Total.Loci=n(),Mean.B=mean(BF), Mean.G=mean(GL), 
            Median.B=median(BF), Median.G=median(GL), 
            Std.B=sd(BF), Std.G=sd(GL),
            Max.B=max(BF), Max.G=max(GL),
            Min.B=min(BF), Min.G=min(GL),
            Out.Sum.B=sum(z.b.outlier), Out.Sum.G=sum(z.g.outlier),
            Out.Prop.B=sum(z.b.outlier)/n(), Out.Prop.G=sum(z.g.outlier)/n(),
            Mad.B=mad(BF), Mad.G=mad(GL), 
            IQR.B=IQR(BF), IQR.G=IQR(GL),
            Skew.B = skewness(BF), Skew.G = skewness(GL), 
            Kurt.B = kurtosis(BF),Kurt.G = kurtosis(GL)) 
#%>% mutate_if(is.numeric,round,20)
# IQR = Q3-Q1 how spread out the middle of our data is

#Outliers based on IQR-- NON STANDARD-----------------------------------------------------------------------------------------------------------------------------------------------------
# Had to alter from <= to < because IQR for GL is 0 in some cases
# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r
# https://www.r-bloggers.com/combined-outlier-detection-with-dplyr-and-ruler/
require("dplyr")

outlier.counter.n <- function(h){
  datf <- all[all$variable==h,]
  d.b <- datf$BF
  d.g <- datf$GL
  b <- IQR.outliers(datf,d.b)
  g <- IQR.outliers(datf,d.g)
  return(c(b,g))
}

# Prints number of outliers according to IQR, bf then gl 
# AIvSA, AIvSI, SAvSI, TvS
stats.all.n$Out.Iqr.B <- c(outlier.counter.n('AIvSA')[2][[1]][2],
                         outlier.counter.n('AIvSI')[2][[1]][2],
                         outlier.counter.n('SAvSI')[2][[1]][2],
                         outlier.counter.n('TvS')[2][[1]][2])

stats.all.n$Out.Iqr.G <- c(outlier.counter.n('AIvSA')[4][[1]][2],
                         outlier.counter.n('AIvSI')[4][[1]][2],
                         outlier.counter.n('SAvSI')[4][[1]][2],
                         outlier.counter.n('TvS')[4][[1]][2])
# Get proportion of loci 
stats.all.n$Out.Iqr.Prop.B <- stats.all.n$Out.Iqr.B/stats.all.n$Total.Loci
stats.all.n$Out.Iqr.Prop.G <- stats.all.n$Out.Iqr.G/stats.all.n$Total.Loci

stats.all.n$variable =paste(stats.all.n$variable, ".og", sep='')

sa <- bind_rows(stats.all,stats.all.n)

# Export to CSV file in folder with "all stats" 

# output excel sheet
library(reshape)
long.stats <- melt(as.data.frame(sa),id=c("variable"))

library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_standardized_distribution.xlsx",sep='')

write.xlsx2(as.data.frame(sa), file=fname, sheetName="ALL", row.names=FALSE)
write.xlsx2(as.data.frame(long.stats), file=fname, sheetName="Long", append=TRUE, row.names=FALSE)

