rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library(praise)
library(svglite)
praise()

# Burbrink ------------------------------------------------------------------------------------------------------------------------------------------------------

dataset <- "Burbrink"

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

# Reshape 
B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

# Smoosh together
a <- bind_rows(B,G) 


# Z score
z.BF <- mLBF %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

z.B <- melt(z.BF,id=c("Locus","supportType")) %>% dplyr::rename(z.value = value)

z.GL <- mLGL %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

z.G <- melt(z.GL,id=c("Locus","supportType")) %>% dplyr::rename(z.value = value)

# Smoosh together
z.a <- bind_rows(z.B,z.G)

# Smoosh together more
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
z.all.e <- left_join(z.G,z.B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
all.ex <- left_join(all.e,z.all.e,by=c("Locus","variable", "supportType.b","supportType.g")) 
all.exp <- all.ex %>% dplyr::rename(BF = value.b,GL = value.g,z.BF = z.value.b,z.GL = z.value.g) %>% 
  select(-c("supportType.b","supportType.g"))

a.z <- left_join(a, z.a, by=c("Locus","variable", "supportType"))

burbrink <- mutate(all.exp, dataSet = case_when(Locus != 'a' ~ dataset))
burbrink.e <- mutate(a.z, dataSet = case_when(Locus != 'a' ~ dataset)) 
rm(all.ex,all.exp,B,G,mLBF,mLGL,a,z.B,z.G,z.BF,z.GL,z.a,dataset, all.e, z.all.e, a.z)

# Reeder ------------------------------------------------------------------------------------------------------------------------------------------------------


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

# Reshape 
B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

# Smoosh together
a <- bind_rows(B,G) 


# Z score
z.BF <- mLBF %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

z.B <- melt(z.BF,id=c("Locus","supportType")) %>% dplyr::rename(z.value = value)

z.GL <- mLGL %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

z.G <- melt(z.GL,id=c("Locus","supportType")) %>% dplyr::rename(z.value = value)

# Smoosh together
z.a <- bind_rows(z.B,z.G)

# Smoosh together more
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
z.all.e <- left_join(z.G,z.B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
all.ex <- left_join(all.e,z.all.e,by=c("Locus","variable", "supportType.b","supportType.g")) 
all.exp <- all.ex %>% dplyr::rename(BF = value.b,GL = value.g,z.BF = z.value.b,z.GL = z.value.g) %>% 
  select(-c("supportType.b","supportType.g"))

a.z <- left_join(a, z.a, by=c("Locus","variable", "supportType"))

reeder <- mutate(all.exp, dataSet = case_when(Locus != 'a' ~ dataset))
reeder.e <- mutate(a.z, dataSet = case_when(Locus != 'a' ~ dataset)) 

rm(all.ex,all.exp,B,G,mLBF,mLGL,a,z.B,z.G,z.BF,z.GL,z.a,dataset, all.e, z.all.e, a.z)

# Singhal ------------------------------------------------------------------------------------------------------------------------------------------------------


dataset <- "Singhal"

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

# Reshape 
B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

# Smoosh together
a <- bind_rows(B,G) 


# Z score
z.BF <- mLBF %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

z.B <- melt(z.BF,id=c("Locus","supportType")) %>% dplyr::rename(z.value = value)

z.GL <- mLGL %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

z.G <- melt(z.GL,id=c("Locus","supportType")) %>% dplyr::rename(z.value = value)

# Smoosh together
z.a <- bind_rows(z.B,z.G)


# Smoosh together more
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
z.all.e <- left_join(z.G,z.B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
all.ex <- left_join(all.e,z.all.e,by=c("Locus","variable", "supportType.b","supportType.g")) 
all.exp <- all.ex %>% dplyr::rename(BF = value.b,GL = value.g,z.BF = z.value.b,z.GL = z.value.g) %>% 
  select(-c("supportType.b","supportType.g"))

a.z <- left_join(a, z.a, by=c("Locus","variable", "supportType"))

singhal <- mutate(all.exp, dataSet = case_when(Locus != 'a' ~ dataset))
singhal.e <- mutate(a.z, dataSet = case_when(Locus != 'a' ~ dataset)) 

rm(all.ex,all.exp,B,G,mLBF,mLGL,a,z.B,z.G,z.BF,z.GL,z.a,dataset, all.e, z.all.e, a.z)

# Streicher ------------------------------------------------------------------------------------------------------------------------------------------------------


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

# Reshape 
B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

# Smoosh together
a <- bind_rows(B,G) 


# Z score
z.BF <- mLBF %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

z.B <- melt(z.BF,id=c("Locus","supportType")) %>% dplyr::rename(z.value = value)

z.GL <- mLGL %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

z.G <- melt(z.GL,id=c("Locus","supportType")) %>% dplyr::rename(z.value = value)

# Smoosh together
z.a <- bind_rows(z.B,z.G)


# Smoosh together more
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
z.all.e <- left_join(z.G,z.B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
all.ex <- left_join(all.e,z.all.e,by=c("Locus","variable", "supportType.b","supportType.g")) 
all.exp <- all.ex %>% dplyr::rename(BF = value.b,GL = value.g,z.BF = z.value.b,z.GL = z.value.g) %>% 
  select(-c("supportType.b","supportType.g"))

a.z <- left_join(a, z.a, by=c("Locus","variable", "supportType"))

streicher <- mutate(all.exp, dataSet = case_when(Locus != 'a' ~ dataset))
streicher.e <- mutate(a.z, dataSet = case_when(Locus != 'a' ~ dataset)) 

rm(all.ex,all.exp,B,G,mLBF,mLGL,a,z.B,z.G,z.BF,z.GL,z.a,dataset, all.e, z.all.e, a.z)

# Standardize Data ------------------------------------------------------------------------------------------------------------------------------------------------------
all.e <- bind_rows(burbrink.e,reeder.e,singhal.e,streicher.e)
all.z <- all.e %>% group_by(dataSet, variable, supportType) %>% mutate_at(c("value"), ~ scale(.,center=TRUE, scale=TRUE))

# Merge Data ------------------------------------------------------------------------------------------------------------------------------------------------------

all <- bind_rows(burbrink,reeder,singhal,streicher)
all <- all %>% 
  mutate(diff = case_when(between(GL,-2,2) & BF > 5 ~ "BF strong dGLS neutral",
                          between(GL,-2,2) & BF < -5 ~ "BF strong dGLS neutral",
                          between(BF,-5,5) & GL > 2 ~ "dGLS strong BF neutral",
                          between(BF,-5,5) & GL < -2 ~ "dGLS strong BF neutral",
                          BF < -5 & GL > 2 ~ "opposite strong support",
                          BF > 5 & GL < -2 ~ "opposite strong support",
                          BF > 5 & GL > 2 ~ "agree",
                          BF < -5 & GL < -2 ~ "agree",
                          between(BF,-5,5) & between(GL,-2,2) ~ "agree"))


all$dataSet <- factor(all$dataSet,levels = c("Reeder","Burbrink","Streicher","Singhal"))
all <- arrange(all, desc(dataSet))

# Subset of Data ------------------------------------------------------------------------------------------------------------------------------------------------------

# group by set numbers - dont need to group by dataset 
all.2x <- all %>% group_by(variable) %>% 
  filter(between(BF,-15,15) | 
           between(GL,-6,6)) # "significance levels (5,2) 5 + 5x 2
all.1x <- all %>% group_by(variable) %>% 
  filter(between(BF,-10,10) | 
           between(GL,-4,4)) # "significance levels (5,2) 5 + 5x 1
all.0x <- all %>% group_by(variable) %>% 
  filter(between(BF,-5,5) | 
           between(GL,-2,2)) # "significance levels (5,2) 5 + 5x 1


# get percentiles - group by dataset to get relative quantiles for each dataset along with each hypoth comparison

all.50l <- all %>% group_by(variable,dataSet) %>% 
  filter(between(GL,quantile(all$GL,0),quantile(all$GL, 0.5)) |
           between(BF,quantile(all$BF,0),quantile(all$BF, 0.5)))

all.50h <- all %>% group_by(variable,dataSet) %>% 
  filter(between(GL,quantile(all$GL,0.5),quantile(all$GL, 1)) |
           between(BF,quantile(all$BF,0.5),quantile(all$BF, 1)))

all.50m <- all %>% group_by(variable,dataSet) %>% 
  filter(between(GL,quantile(all$GL,0.25),quantile(all$GL, 0.75)) |
           between(BF,quantile(all$BF,0.25),quantile(all$BF, 0.75)))

# all bf within and all gl within. agnostic to other value. gl in OR bf in 
all.60m <- all %>% group_by(variable,dataSet) %>% 
  filter(between(GL,quantile(all$GL,0.2),quantile(all$GL, 0.8)) | 
           between(BF,quantile(all$BF,0.2),quantile(all$BF, 0.8)))

all.80m <- all %>% group_by(variable,dataSet) %>% 
  filter(between(GL,quantile(all$GL,0.1),quantile(all$GL, 0.9)) | 
           between(BF,quantile(all$BF,0.1),quantile(all$BF, 0.9)))

# both within - gl AND bf within  
all.60mb <- all %>% group_by(variable,dataSet) %>% 
  filter(between(GL,quantile(all$GL,0.2),quantile(all$GL, 0.8)) , 
           between(BF,quantile(all$BF,0.2),quantile(all$BF, 0.8)))

all.80mb <- all %>% group_by(variable,dataSet) %>% 
  filter(between(GL,quantile(all$GL,0.1),quantile(all$GL, 0.9)) , 
           between(BF,quantile(all$BF,0.1),quantile(all$BF, 0.9)))



# Save Data ------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
save.image(file=paste("Calcs_Smoosh.RData",sep=""))

# # Disagree Data ------------------------------------------------------------------------------------------------------------------------------------------------------
# # both sig
# # both same direction
# # both opposite direction
# # both sig and opposite 
# 
# g <- 2
# b <- 5
# #d <- 0.5
# # grouping when both are 0
# # when 1 is zero, counting as "same"
# s <- all %>% dplyr::rename(Hypothesis=variable) %>% 
#   mutate(sig = case_when(BF==0 & GL==0 ~ 'both.zero.zero',
#                          BF==0 & GL!=0 ~ 'BF.zero.zero',
#                          BF!=0 & GL==0 ~ 'GL.zero.zero',
#                          BF >= b & GL >= g ~ 'both.sig.same',
#                          BF <= -b & GL <= -g ~ 'both.sig.same',
#                          BF <= -b & GL >= g ~ 'both.sig.opposite',
#                          BF >= b & GL <= -g ~ 'both.sig.opposite',
#                          BF > 0 & GL > 0 ~ 'both.any.same',
#                          BF < 0 & GL < 0 ~ 'both.any.same',
#                          BF > 0 & GL < 0 ~ 'both.any.opposite',
#                          BF < 0 & GL > 0 ~ 'both.any.opposite'
#   )) 
# 
# s %>% group_by(sig) %>% summarise(n=n())
# 
# 
# %>% 
#   separate(sig,c("a","b","c"),remove=FALSE) %>% 
#   mutate(d=case_when(c=='same' & b == 'sigB' ~ 'same, one sig',
#                      c=='same' & b == 'sigG' ~ 'same, one sig',
#                      c=='same' & b == 'neu' ~ 'same, both neu',
#                      c=='same' & b == 'sig' ~ 'same, both sig',
#                      c=='opposite' ~ 'opposite',
#                      c=='zero' ~ 'zero')) %>%
#   mutate(Hypothesis=factor(Hypothesis, levels=c("AIvSA","AIvSI","SAvSI","TvS"))) %>%
#   mutate(c=factor(c,levels = c("zero","opposite","same"))) %>%
#   mutate(h=case_when(Hypothesis=='TvS' ~ 'TvS',
#                      Hypothesis!='TvS' ~ 'Tox')) %>% 
#   mutate(thresh=case_when(a == 'neuB' ~ 'conflict',
#                           a == 'neuG' ~ 'conflict',
#                           sig == 'both.sig.opposite' ~ 'conflict',
#                           sig == 'both.neu.opposite' ~ 'agree',
#                           sig == 'both.neu.same' ~ 'agree',
#                           sig == 'both.sig.same' ~ 'agree',
#                           sig == 'both.zero.zero' ~ 'agree'))
# 
# 
# # Count loci that have BF and dGLS that are in the same direction of support and those in opposite, and those that are zero
# support.direction <- s %>% dplyr::group_by(c,Hypothesis,.drop=FALSE) %>% dplyr::summarise(n.loci=n(),.groups = 'keep') %>% 
#   mutate(percent.loci= round((n.loci/loci)*100,digits = 2)) %>% 
#   ungroup() %>% group_by(Hypothesis) %>%  mutate(a=sum(n.loci),b=sum(percent.loci)) 
# 
# 
# # Count loci that disagree in terms of thresholds 
# t.hold <- s %>% dplyr::group_by(thresh,Hypothesis,.drop=FALSE) %>% dplyr::summarise(n.loci=n(),.groups = 'keep') %>% 
#   mutate(percent.loci= round((n.loci/loci)*100,digits = 2)) %>% 
#   ungroup() %>% group_by(Hypothesis) %>%  mutate(a=sum(n.loci),b=sum(percent.loci)) %>%
#   mutate(thresh=factor(thresh,levels = c("conflict","agree")))
# 
# 
# # Tally all loci where significance thresholds conflict - both sig and opposite. one sig and opp, one sig and same
# t.hold.details <- s  %>% filter(!sig %in% c('both.neu.same','both.neu.opposite','both.zero.zero','both.sig.same')) %>%
#   separate(sig,c("a","b","c"),remove=FALSE) %>%
#   unite("ab",a:b, remove=FALSE) %>%
#   dplyr::group_by(ab,Hypothesis,.drop=FALSE) %>% dplyr::summarise(n=n(),.groups = 'keep') %>% 
#   ungroup() %>% group_by(Hypothesis) %>%
#   dplyr::mutate(percent.loci= round((n/loci)*100,digits = 2)) %>% 
#   droplevels() %>% mutate(thresh=factor(ab,levels = c("both_sig","neuG_sigB","neuB_sigG")))
# 
# 
# 
# all.d <- all %>% 
#   mutate(diff = case_when(between(GL,-2,2) & BF > 5 ~ "BF strong dGLS neutral",
#                           between(GL,-2,2) & BF < -5 ~ "BF strong dGLS neutral",
#                           between(BF,-5,5) & GL > 2 ~ "dGLS strong BF neutral",
#                           between(BF,-5,5) & GL < -2 ~ "dGLS strong BF neutral",
#                           BF < -5 & GL > 2 ~ "opposite strong support",
#                           BF > 5 & GL < -2 ~ "opposite strong support",
#                           BF > 5 & GL > 2 ~ "agree",
#                           BF < -5 & GL < -2 ~ "agree",
#                           between(BF,-5,5) & between(GL,-2,2) ~ "agree"))

