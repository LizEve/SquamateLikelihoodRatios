rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)


##################################################################################################
## set dataset 
##################################################################################################



dataset <- "Reeder"

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
all <- bind_rows(B,G) 


all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
all.exp <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% select(-c("supportType.b","supportType.g"))


all.60m <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) | 
           between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80m <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) | 
           between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# both within - gl AND bf within  
all.60mb <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) , 
         between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80mb <- all.exp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) , 
         between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

##################################################################################################
## equations
##################################################################################################


#y-gl
#x-bf
lr.df<- data.frame(hypothesis=character(),
                   subset=character(),
                   a=numeric(), 
                   b=numeric(), 
                   r2=numeric(),
                   nloci=numeric(),
                   source=character())
dsets <- list(all=all.exp,all.80m=all.80m,all.60m=all.60m,all.80mb=all.80mb,all.60mb=all.60mb) 

hypos <- c("TvS","AIvSA","AIvSI","SAvSI")

for (x in seq(1,length(dsets))){
  dd <- names(dsets[x])
  d <- as.data.frame(dsets[x])
  v <- paste(dd,"variable",sep=".")
  print(d)
  print(v)
  for (h in hypos){
    print(dd)
    print(h)
    df <- d %>% filter(d[[v]] == h)
    print(df)
    b <- paste(dd,"BF",sep=".")
    g <- paste(dd,"GL",sep=".")
    m = lm(df[[b]] ~ df[[g]])
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


# output excel sheet

library(xlsx)


fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/Tables/LinearRegression_",dataset,"_mb_2021.xlsx",sep='')

write.xlsx2(as.data.frame(lr.df), file=fname, sheetName=dataset, row.names=FALSE,append=TRUE)


