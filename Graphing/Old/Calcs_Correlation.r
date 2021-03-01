rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)

# set dataset 

dataset <- "Burbrink"

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
#mLBF[4] <- mLBF[4]/2

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


keep <- c(1,seq(14,22))

B <- bind_rows(reformatdf(mLBF,"BF","AIvSA",1,keep),
               reformatdf(mLBF,"BF","AIvSI",1,keep),
               reformatdf(mLBF,"BF","SAvSI",1,keep),
               reformatdf(mLBF,"BF","TvS",1,keep))

G <- bind_rows(reformatdf(mLGL,"GL","AIvSA",1,keep),
               reformatdf(mLGL,"GL","AIvSI",1,keep),
               reformatdf(mLGL,"GL","SAvSI",1,keep),
               reformatdf(mLGL,"GL","TvS",1,keep))

# for SinghalOG 
#keep <- c(1,seq(5,13)) 
#B <- reformatdf(mLBF,"BF","TvS",1,keep)
#G <- reformatdf(mLGL,"GL","TvS",1,keep)


metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")
all <- left_join(B, G, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))

##################################################################################################
## Test correlation across all hypotheses 
##################################################################################################
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# https://rpubs.com/aaronsc32/spearman-rank-correlation
# Run from current line to document end 	Ctrl+Alt+E

rm(null,null_rs,new,sig2tp05,emp_rs,a,h)
h <- 'TvS'
#h <- 'AIvSA'
#h <- 'AIvSI'
#h <- 'SAvSI'
a <- all[all$Hypothesis == h,]

# start with all the meta data 

null <- a %>% select(c(Locus,Hypothesis,metaData,BF,GL))

# Permute BF and GL in 1,000 seperate columns

for (n in 1:1000) {
  newcolname1 <- paste('BF',n,sep='') 
  null[, newcolname1] <- sample(a$BF, replace = FALSE)
  newcolname2 <- paste('GL',n,sep='') 
  null[, newcolname2] <- sample(a$GL, replace = FALSE)
}

# Create new dataset for spearman rank coefficients

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

save(all,null,null_rs, file=paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_CorrNullDist_",h,".RData",sep=''))
#save(null_rs, file=paste("/Users/ChatNoir/Projects/Squam/scripts/Graphing/",dataset,"_CorrNullDist_",h,"_rs.RData",sep=''))

#load("/Users/ChatNoir/Projects/Squam/scripts/Graphing/CorrelationNullDist_rs.RData")

# Add row of emp scores 

assign('emp_rs', setNames(data.frame('emp'), "replicate"))
empRow <- 1

# First bf and gl 

# Grab column names, calculate value for BF GL there is only one value for both columns

i <- "BF"
cc1 <- cor(a[,i],a[,'GL'],method="spearman") 
empcol1 <- paste(i,'.g.emp',sep='')
emp_rs[empRow,empcol1] <- cc1
i <- "GL"
empcol2 <- paste(i,'.b.emp',sep='')
emp_rs[empRow,empcol2] <- cc1

# Add all metadata emp for BF and GL 

metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")

for (i in metaData) {
  cc1 <- cor(a[,i],a[,'GL'],method="spearman") 
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

# Save null and emp rs as backup 
save(null_rs, emp_rs, sig2tp05,  file=paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",h,"_rsvals.RData",sep=''))

# Save each hypothesis to output in xlxs sheet for each dataset 

assign(h, sig2tp05)


#ctrl shift c - to uncomment. then run at end
#https://stackoverflow.com/questions/48799959/r-reduce-with-merge-and-more-than-2-suffixes-or-how-to-merge-multiple-datafr
#Merge all datasets

list.df <- list(TvS,AIvSA,AIvSI,SAvSI)
sfx <- c("_TvS", "_AIvSA", "_AIvSI", "_SAvSI")
res <- list.df[[1]]
for(i in head(seq_along(list.df), -1)) {

  res <- merge(res, list.df[[i+1]], all = TRUE,
               suffixes = sfx[i:(i+1)], by = c("metaData"))
}

# output excel sheet

library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_correlation.xlsx",sep='')

write.xlsx2(as.data.frame(res), file=fname, sheetName="ALL", row.names=FALSE)
write.xlsx2(as.data.frame(TvS), file=fname, sheetName="TvS", append=TRUE, row.names=FALSE)
write.xlsx2(as.data.frame(AIvSA), file=fname, sheetName="AIvSA", append=TRUE, row.names=FALSE)
write.xlsx2(as.data.frame(AIvSI), file=fname, sheetName="AIvSI", append=TRUE, row.names=FALSE)
write.xlsx2(as.data.frame(SAvSI), file=fname, sheetName="SAvSI", append=TRUE, row.names=FALSE)

