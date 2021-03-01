rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)


# set dataset 

dataset <- "Streicher"

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
all <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% select(-c("supportType.b","supportType.g"))

# Subset of Data ------------------------------------------------------------------------------------------------------------------------------------------------------
# Originally ran with threshold of dGLS being 0.5, 2x = 1, 4x=2
# trying again with 2, 2x=4, 4x=8

all.sx4 <- all %>% group_by(variable) %>% 
  filter(between(BF,-20,20)) %>% 
  filter(between(GL,-8,8)) # "significance levels (5,0.5) x 4
all.sx2 <- all %>% group_by(variable) %>% 
  filter(between(BF,-10,10)) %>% 
  filter(between(GL,-4,4)) # "significance levels (5,0.5) x 2
all.mid <- all %>% group_by(variable) %>% 
  filter(between(BF,quantile(all$BF, 0.25),quantile(all$BF, 0.75))) %>% 
  filter(between(GL,quantile(all$GL, 0.25),quantile(all$GL, 0.75))) # Middle 50% of distribution

# all.mid <- all %>% group_by(variable) %>% 
#   filter(between(BF,quantile(all$BF, 0.1),quantile(all$BF, 0.9))) %>% 
#   filter(between(GL,quantile(all$GL, 0.1),quantile(all$GL, 0.9))) # Middle 50% of distribution



#y-gl
#x-bf
lr.df<- data.frame(hypothesis=character(),
                   subset=character(),
                   a=numeric(), 
                   b=numeric(), 
                   r2=numeric(),
                   nloci=numeric(),
                   source=character())
dsets <- list(all=all,all.sx4=all.sx4,all.sx2=all.sx2,all.mid=all.mid) 
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

rm(d,dd,v,g,x,m,l)

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

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_linearRegression2.xlsx",sep='')

write.xlsx2(as.data.frame(lr.df), file=fname, sheetName="LinearRegression2", row.names=FALSE)
