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

#a.z <- left_join(a, z.a, by=c("Locus","variable", "supportType"))

# add column with dataset name 
temp <- mutate(all.exp, dataSet = case_when(Locus != 'a' ~ dataset))


# all bf within and all gl within. agnostic to other value. gl in OR bf in 
all.60m <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) | 
           between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80m <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) | 
           between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# both within - gl AND bf within  
all.60mb <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) , 
         between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80mb <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) , 
         between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# calculate rank - removing everything except what im graphing for now
comp <- temp %>% dplyr::group_by(variable) %>% 
  #mutate(per.g = ntile(desc(GL),100)) %>% 
  mutate(drank.g = dense_rank(desc(GL))) %>%
  #mutate(rank.g = rank(desc(GL))) %>%
  #mutate(per.b = ntile(desc(BF),100)) %>% 
  mutate(drank.b = dense_rank(desc(BF))) %>%
  #mutate(rank.b = rank(desc(BF))) %>%
  #mutate(per.diff = abs(per.g-per.b)) %>%
  #mutate(drank.diff = abs(drank.g - drank.b)) %>%
  mutate(drank.diff.n = abs((drank.g/max(drank.g))-(drank.b/max(drank.b)))*100) #%>%
  #mutate(rank.diff = abs(rank.g - rank.b)) %>%
  #mutate(rank.diff.n = abs((rank.g/max(rank.g))-(rank.b/max(rank.b)))*100)


# reassign variables with dataset specific names 
burbrink <- comp
bur.60m <- all.60m
bur.80m <- all.80m
bur.60mb <- all.60mb
bur.80mb <- all.80mb
#burbrink.e <- mutate(a.z, dataSet = case_when(Locus != 'a' ~ dataset)) 
rm(comp,a,all.60m,all.60mb,all.80m, all.80mb,all.e,all.exp,B,G,mLBF, mLGL,temp,z.a, z.all.e,z.B,z.BF,z.G,z.GL,all.ex)

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

# add column with dataset name 
temp <- mutate(all.exp, dataSet = case_when(Locus != 'a' ~ dataset))


# all bf within and all gl within. agnostic to other value. gl in OR bf in 
all.60m <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) | 
           between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80m <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) | 
           between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# both within - gl AND bf within  
all.60mb <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) , 
         between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80mb <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) , 
         between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# calculate rank - removing everything except what im graphing for now
comp <- temp %>% dplyr::group_by(variable) %>% 
  mutate(drank.g = dense_rank(desc(GL))) %>%
  mutate(drank.b = dense_rank(desc(BF))) %>%
  mutate(drank.diff.n = abs((drank.g/max(drank.g))-(drank.b/max(drank.b)))*100) 

# reassign variables with dataset specific names 
reeder <- comp
ree.60m <- all.60m
ree.80m <- all.80m
ree.60mb <- all.60mb
ree.80mb <- all.80mb
rm(comp,a,all.60m,all.60mb,all.80m, all.80mb,all.e,all.exp,B,G,mLBF, mLGL,temp,z.a, z.all.e,z.B,z.BF,z.G,z.GL,all.ex)

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

# add column with dataset name 
temp <- mutate(all.exp, dataSet = case_when(Locus != 'a' ~ dataset))


# all bf within and all gl within. agnostic to other value. gl in OR bf in 
all.60m <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) | 
           between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80m <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) | 
           between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# both within - gl AND bf within  
all.60mb <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) , 
         between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80mb <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) , 
         between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# calculate rank - removing everything except what im graphing for now
comp <- temp %>% dplyr::group_by(variable) %>% 
  mutate(drank.g = dense_rank(desc(GL))) %>%
  mutate(drank.b = dense_rank(desc(BF))) %>%
  mutate(drank.diff.n = abs((drank.g/max(drank.g))-(drank.b/max(drank.b)))*100) 

# reassign variables with dataset specific names 
singhal <- comp
sin.60m <- all.60m
sin.80m <- all.80m
sin.60mb <- all.60mb
sin.80mb <- all.80mb
rm(comp,a,all.60m,all.60mb,all.80m, all.80mb,all.e,all.exp,B,G,mLBF, mLGL,temp,z.a, z.all.e,z.B,z.BF,z.G,z.GL,all.ex)

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

# add column with dataset name 
temp <- mutate(all.exp, dataSet = case_when(Locus != 'a' ~ dataset))


# all bf within and all gl within. agnostic to other value. gl in OR bf in 
all.60m <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) | 
           between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80m <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) | 
           between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# both within - gl AND bf within  
all.60mb <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) , 
         between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80mb <- temp %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) , 
         between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# calculate rank - removing everything except what im graphing for now
comp <- temp %>% dplyr::group_by(variable) %>% 
  mutate(drank.g = dense_rank(desc(GL))) %>%
  mutate(drank.b = dense_rank(desc(BF))) %>%
  mutate(drank.diff.n = abs((drank.g/max(drank.g))-(drank.b/max(drank.b)))*100) 

# reassign variables with dataset specific names 
streicher <- comp
str.60m <- all.60m
str.80m <- all.80m
str.60mb <- all.60mb
str.80mb <- all.80mb
rm(comp,a,all.60m,all.60mb,all.80m, all.80mb,all.e,all.exp,B,G,mLBF, mLGL,temp,z.a, z.all.e,z.B,z.BF,z.G,z.GL,all.ex)

# Merge Data ------------------------------------------------------------------------------------------------------------------------------------------------------

all <- bind_rows(burbrink,reeder,singhal,streicher)
all.60m <- bind_rows(bur.60m,ree.60m,sin.60m,str.60m)
all.80m <- bind_rows(bur.80m,ree.80m,sin.80m,str.80m)
all.60mb <- bind_rows(bur.60mb,ree.60mb,sin.60mb,str.60mb)
all.80mb <- bind_rows(bur.80mb,ree.80mb,sin.80mb,str.80mb)

all$dataSet <- factor(all$dataSet,levels = c("Reeder","Burbrink","Streicher","Singhal"))
all.60m$dataSet <- factor(all.60m$dataSet,levels = c("Reeder","Burbrink","Streicher","Singhal"))
all.80m$dataSet <- factor(all.80m$dataSet,levels = c("Reeder","Burbrink","Streicher","Singhal"))
all.60mb$dataSet <- factor(all.60mb$dataSet,levels = c("Reeder","Burbrink","Streicher","Singhal"))
all.80mb$dataSet <- factor(all.80mb$dataSet,levels = c("Reeder","Burbrink","Streicher","Singhal"))
all <- arrange(all, desc(dataSet))


# Save Data ------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
save.image(file=paste("Calcs_Smoosh2021April.RData",sep=""))
