rm(list=ls())
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

B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))


G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

# Smoosh together again
a <- bind_rows(B,G) 
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
all.exp <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% select(-c("supportType.b","supportType.g"))
all <- all.exp
# Subset of Data ------------------------------------------------------------------------------------------------------------------------------------------------------


all.2x <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-15,15) | 
           between(GL,-6,6)) # "significance levels (5,2) 5 + 5x 2
all.1x <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-10,10) | 
           between(GL,-4,4)) # "significance levels (5,2) 5 + 5x 1
all.0x <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-5,5) | 
           between(GL,-2,2)) # "significance levels (5,2) 5 + 5x 1


# get percentiles 

all.50l <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0),quantile(all.exp$GL, 0.5)) |
           between(BF,quantile(all.exp$BF,0),quantile(all.exp$BF, 0.5)))

all.50h <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0.5),quantile(all.exp$GL, 1)) |
           between(BF,quantile(all.exp$BF,0.5),quantile(all.exp$BF, 1)))

all.50m <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0.25),quantile(all.exp$GL, 0.75)) |
           between(BF,quantile(all.exp$BF,0.25),quantile(all.exp$BF, 0.75)))

all.60m <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0.2),quantile(all.exp$GL, 0.8)) | 
           between(BF,quantile(all.exp$BF,0.2),quantile(all.exp$BF, 0.8)))

all.80m <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.exp$GL,0.1),quantile(all.exp$GL, 0.9)) | 
           between(BF,quantile(all.exp$BF,0.1),quantile(all.exp$BF, 0.9)))




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
              all.80m=all.80m,all.60m=all.60m,all.50m=all.50m,
              all.50l=all.50l,all.50h=all.50h) 
#dsets <- list(all=all,all.sx4=all.sx4,all.sx2=all.sx2) 

hypos <- c("TvS","AIvSA","AIvSI","SAvSI")
#hypos <- c("AIvSA","AIvSI","SAvSI") # reeder

# # TvS only appears in all dataset for reeder
# df <- all %>% filter(variable == "TvS")
# m = lm(df$GL ~ df$BF)
# l <- list(hypothesis="TvS", 
#          subset="all", 
#          a = format(unname(coef(m)[1]), digits = 2),
#          b = format(unname(coef(m)[2]), digits = 2), 
#          r2 = format(summary(m)$r.squared, digits = 3),
#          nloci = length(df$GL),
#          source = dataset)
# lr.df <- rbind(lr.df,data.frame(l))


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
              nloci = length(df[[g]]),
              source = dataset)
    #print(data.frame(l))
    lr.df <- rbind(lr.df,data.frame(l))
  }
}


library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_linearRegression_Sept16inclusive.xlsx",sep='')

write.xlsx2(as.data.frame(lr.df), file=fname, sheetName="LinearRegression2", row.names=FALSE)

