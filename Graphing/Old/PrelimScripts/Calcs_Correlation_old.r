rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)

# set dataset 

dataset <- "SinghalOG"

# Colors 
# https://www.w3schools.com/colors/colors_picker.asp
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS"))
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10
mLBF[7:10] <- mLBF[7:10]/2
# for singhal OG 
mLBF[4] <- mLBF[4]/2

# Check number of NA 
sum(is.na(mLBF$MEAN_COL_SCORE))

# Transform datasets - stack two comparisons, and add column for bf and gl 
# Each dataset one hypothesis 

reformatdf <- function(df,h,h1,s1,keepCols){
  # keepCols is a vector of the locus and support type column indices
  # Grab column number  
  c1 <- match(h1,names(df))
  # Grab and rename column, add hypothesis column
  a.df <- df[,append(keepCols,c1)] 
  names(a.df)[names(a.df) == h1] <- h
  a.df$Hypothesis <- rep(h1,length(a.df[,1]))
  c <- match(h,names(a.df))
  # Adjust direction of support if needed. ie if column is AIvSA, but you want to know SAvAI
  a.df[,c] <- a.df[,c]*(s1)
  return(a.df)
}

# reformat 
keep <- c(1,seq(14,22))
AI.g <- bind_rows(reformatdf(mLGL,"AI","AIvSA",1,keep),
                  reformatdf(mLGL,"AI","AIvSI",1,keep))

SA.g <- bind_rows(reformatdf(mLGL,"SA","AIvSA",-1,keep),
                  reformatdf(mLGL,"SA","SAvSI",1,keep))

SI.g <- bind_rows(reformatdf(mLGL,"SI","AIvSI",-1,keep),
                  reformatdf(mLGL,"SI","SAvSI",-1,keep))

AI.b <- bind_rows(reformatdf(mLBF,"AI","AIvSA",1,keep),
                  reformatdf(mLBF,"AI","AIvSI",1,keep))

SA.b <- bind_rows(reformatdf(mLBF,"SA","AIvSA",-1,keep),
                  reformatdf(mLBF,"SA","SAvSI",1,keep))

SI.b <- bind_rows(reformatdf(mLBF,"SI","AIvSI",-1,keep),
                  reformatdf(mLBF,"SI","SAvSI",-1,keep))

keep <- c(1,seq(5,13)) # for SinghalOG 
TS.g <- reformatdf(mLGL,"TS","TvS",-1,keep)
TS.b <- reformatdf(mLBF,"TS","TvS",-1,keep)

# MERGE - yea i could have merged first, but the function was already written and i didnt want to mess with it. 

AI.a <- left_join(AI.g, AI.b, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))
SA.a <- left_join(SA.g, SA.b, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))
SI.a <- left_join(SI.g, SI.b, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))
TS.a <- left_join(TS.g, TS.b, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))

sapply(AI.a,typeof)

# Smash all together, then smash bf and gl together each in 1 column
all <- bind_rows(bind_rows(bind_rows(AI.a,SA.a),SI.a),TS.a) %>% 
  mutate(BF=coalesce(AI.b,SA.b,SI.b,TS.b)) %>% 
  mutate(GL=coalesce(AI.g,SA.g,SI.g,TS.g)) %>%
  select(-c(AI.b,SA.b,SI.b,TS.b,AI.g,SA.g,SI.g,TS.g))

TX <- bind_rows(bind_rows(AI.a,SA.a),SI.a) %>% 
  mutate(BF=coalesce(AI.b,SA.b,SI.b)) %>% 
  mutate(GL=coalesce(AI.g,SA.g,SI.g)) %>%
  select(-c(AI.b,SA.b,SI.b,AI.g,SA.g,SI.g))


TS <- TS.a %>% 
  mutate(BF=coalesce(TS.b)) %>% 
  mutate(GL=coalesce(TS.g)) %>%
  select(-c(TS.b,TS.g))


# Check number of NA 
sum(is.na(all$BF))
sum(is.na(all$GL))

rm(AI.g,AI.b,SA.b,SA.g,SI.b,SI.g,TS.b,TS.g)
rm(AI.a,SA.a,SI.a,TS.a)
rm(mLBF,mLGL)
##################################################################################################
## Test correlation across all hypotheses 
##################################################################################################
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# https://rpubs.com/aaronsc32/spearman-rank-correlation
rm(null,null_rs,new,sig2tp05,emp_rs)
# Create null distribution 
all <- TS
v <- "ts"
# start with all the meta data 

null <- all %>% select(c(Locus,Hypothesis,metaData,BF,GL))

# Permute BF and GL in 1,000 seperate columns

for (n in 1:1000) {
  newcolname1 <- paste('BF',n,sep='') 
  null[, newcolname1] <- sample(all$BF, replace = FALSE)
  newcolname2 <- paste('GL',n,sep='') 
  null[, newcolname2] <- sample(all$GL, replace = FALSE)
}

# Create new dataset for spearman rank coefficients (repeat later with kendall?)

assign('null_rs', setNames(data.frame(1:1000), "replicate"))

# for meta data 

for (x in metaData) {
  null_rs[,paste(x,'.b.null',sep='')] <- NA
  null_rs[,paste(x,'.g.null',sep='')] <- NA
}

# then for bf and gl - these should end up the same 

null_rs[,'GL.b.null'] <- NA # bf replicated, gl same
null_rs[,'BF.g.null'] <- NA # gl replicated, bf same

# Add bf gl 

for (n in 1:1000) {
  
  # Grab column names, calculate value 
  colname1 <- paste('BF',n,sep='') 
  cc1 <- cor(null[,'GL'],null[,colname1],method="spearman") 
  nullcol1 <- paste('GL','.b.null',sep='')
  null_rs[n,nullcol1] <- cc1
  colname2 <- paste('GL',n,sep='') 
  cc2 <- cor(null[,'BF'],null[,colname2],method="spearman") 
  nullcol2 <- paste('BF','.g.null',sep='')
  null_rs[n,nullcol2] <- cc2
  
}


for (i in metaData) {
  
  # Calculate rs for all permuted columns in null dataset
  
  for (n in 1:1000) {
    
    # Grab column names, calculate value 
    colname1 <- paste('BF',n,sep='') 
    cc1 <- cor(null[,i],null[,colname1],method="spearman") 
    nullcol1 <- paste(i,'.b.null',sep='')
    null_rs[n,nullcol1] <- cc1
    colname2 <- paste('GL',n,sep='') 
    cc2 <- cor(null[,i],null[,colname2],method="spearman") 
    nullcol2 <- paste(i,'.g.null',sep='')
    null_rs[n,nullcol2] <- cc2
    
  }
}

save(all,null,null_rs, file=paste("/Users/ChatNoir/Projects/Squam/scripts/Graphing/",dataset,"_CorrNullDist_",v,".RData",sep=''))
save(null_rs, file=paste("/Users/ChatNoir/Projects/Squam/scripts/Graphing/",dataset,"_CorrNullDist_",v,"_rs.RData",sep=''))

#load("/Users/ChatNoir/Projects/Squam/scripts/Graphing/CorrelationNullDist_rs.RData")

# Add row of emp scores 

assign('emp_rs', setNames(data.frame('emp'), "replicate"))
empRow <- 1

# First bf and gl 

# Grab column names, calculate value for BF GL there is only one value for both columns

i <- "BF"
cc1 <- cor(all[,i],all[,'GL'],method="spearman") 
empcol1 <- paste(i,'.g.emp',sep='')
emp_rs[empRow,empcol1] <- cc1
i <- "GL"
empcol2 <- paste(i,'.b.emp',sep='')
emp_rs[empRow,empcol2] <- cc1

# Add all metadata emp for BF and GL 

metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")

for (i in metaData) {
  cc1 <- cor(all[,i],all[,'GL'],method="spearman") 
  empcol1 <- paste(i,'.g.emp',sep='')
  emp_rs[empRow,empcol1] <- cc1
  cc2 <- cor(all[,i],all[,'BF'],method="spearman") 
  empcol2 <- paste(i,'.b.emp',sep='')
  emp_rs[empRow,empcol2] <- cc2
}

# Add rows for quantile cut offs 

emp_rs$replicate <- as.character(emp_rs$replicate)
emp_rs[2,1] <- "p0.025"
emp_rs[3,1] <- "p0.975"

# add gl and bf 

emp_rs[2,'GL.b.emp'] <- quantile(null_rs$BF.g.null, c(0.025, 0.975))[[1]]
emp_rs[3,'GL.b.emp'] <- quantile(null_rs$BF.g.null, c(0.025, 0.975))[[2]]
emp_rs[2,'BF.g.emp'] <- quantile(null_rs$GL.b.null, c(0.025, 0.975))[[1]]
emp_rs[3,'BF.g.emp'] <- quantile(null_rs$GL.b.null, c(0.025, 0.975))[[2]]

for (i in metaData) {
  gl.n <- paste(i,'.g.null',sep='')
  gl.e <- paste(i,'.g.emp',sep='')
  g025 <- quantile(null_rs[,gl.n], c(0.025, 0.975))[[1]]
  g975 <- quantile(null_rs[,gl.n], c(0.025, 0.975))[[2]]
  emp_rs[2,gl.e] <- g025
  emp_rs[3,gl.e] <- g975
  bf.n <- paste(i,'.b.null',sep='')
  bf.e <- paste(i,'.b.emp',sep='')
  b025 <- quantile(null_rs[,bf.n], c(0.025, 0.975))[[1]]
  b975 <- quantile(null_rs[,bf.n], c(0.025, 0.975))[[2]]
  emp_rs[2,bf.e] <- b025
  emp_rs[3,bf.e] <- b975
}

# Create two boolean rows for significance 

# Transform dataset from rows to columns 

new <- data.frame(t(emp_rs))
colnames(new) <- as.character(unlist(new[1,]))
new <- new[-1,]
new$emp <- as.numeric(as.character(new$emp))
new$p0.025 <- as.numeric(as.character(new$p0.025))
new$p0.975 <- as.numeric(as.character(new$p0.975))
sig2tp05 <- new %>% rownames_to_column('metaData') %>%
  mutate(sig = case_when(emp < p0.025 ~ "Lower",
                                       emp > p0.975 ~ "Upper",
                                    emp >= p0.025 & emp <= p0.975 ~ "False"))

# Save null and emp rs 
save(null_rs, emp_rs, sig2tp05,  file=paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",v,"_rsvals.RData",sep=''))

# Load and save as csv 
dataset <- "Singhal"
v <- 'all'
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",v,"_rsvals.RData",sep=''))
write.csv(sig2tp05,paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",v,"_rsvals.csv",sep=''))
rm(v,sig2tp05)
v <- 'tx'
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",v,"_rsvals.RData",sep=''))
write.csv(sig2tp05,paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",v,"_rsvals.csv",sep=''))
rm(v,sig2tp05)
v <- 'ts'
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",v,"_rsvals.RData",sep=''))
write.csv(sig2tp05,paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",v,"_rsvals.csv",sep=''))
rm(dataset,v,sig2tp05)
