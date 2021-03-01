rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)

# set dataset 

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

# Check number of NA 
sum(is.na(mLBF$MEAN_COL_SCORE))

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

# Subset data 
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 


all.exp <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% select(-c("supportType.b","supportType.g"))
all.short <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-20,20)) %>% 
  filter(between(GL,-8,8)) # "significance levels (5,0.5) x 4
all.shorter <- all.exp %>% group_by(variable) %>% 
  filter(between(BF,-10,10)) %>% 
  filter(between(GL,-8,8)) # "significance levels (5,0.5) x 2
# Filter based on only one threshold
#all.short.b <- all.exp %>% group_by(variable) %>% 
#  filter(between(BF,-20,20)) # "significance levels (5,0.5) x 4
#all.shorter.b <- all.exp %>% group_by(variable) %>% 
#  filter(between(BF,-10,10))  # "significance levels (5,0.5) x 2
#all.short.g <- all.exp %>% group_by(variable) %>% 
#  filter(between(GL,-2,2)) # "significance levels (5,0.5) x 4
#all.shorter.g <- all.exp %>% group_by(variable) %>% 
#  filter(between(GL,-4,4)) # "significance levels (5,0.5) x 2
#all.mid.b <- all.exp %>% group_by(variable) %>% 
#  filter(between(BF,quantile(all.exp$BF, 0.25),quantile(all.exp$BF, 0.75)))  # Middle 50% of distribution
#all.mid.g <- all.exp %>% group_by(variable) %>% 
#  filter(between(GL,quantile(all.exp$GL, 0.25),quantile(all.exp$GL, 0.75))) # Middle 50% of distribution


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



null <- a 

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

#for (x in metaData) {
#  null_rs[,paste(x,'.b.null',sep='')] <- NA
#  null_rs[,paste(x,'.g.null',sep='')] <- NA
#}

# then for bf and gl - these should end up the same 

null_rs[,'GL.b.null'] <- NA # bf replicated, gl same
null_rs[,'BF.g.null'] <- NA # gl replicated, bf same

# Add bf gl 

for (n in 1:1000) {
  
  # Grab column names, calculate value 
  colname1 <- paste('BF',n,sep='') 
  cc1.s <- cor(null[,'GL'],null[,colname1],method="spearman") 
  null_rs[n,'GL.b.null.spearman'] <- cc1.s
  cc1.p <- cor(null[,'GL'],null[,colname1],method="pearson") 
  null_rs[n,'GL.b.null.pearson'] <- cc1.p
  colname2 <- paste('GL',n,sep='') 
  cc2.s <- cor(null[,'BF'],null[,colname2],method="spearman") 
  null_rs[n,'BF.g.null.spearman'] <- cc2.s
  cc2.p <- cor(null[,'BF'],null[,colname2],method="pearson") 
  null_rs[n,'BF.g.null.pearson'] <- cc2.p
  
}


#save(all,null,null_rs, file=paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_BFGLNullDist_",h,".RData",sep=''))

#load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_BFGLNullDist_",h,".RData",sep=''))

# Add row of emp scores 

assign('emp_rs', setNames(data.frame('emp'), "replicate"))
empRow <- 1

# First bf and gl 

# Grab column names, calculate value for BF GL there is only one value for both columns

i <- "BF"
cc1.s <- cor(a[,i],a[,'GL'],method="spearman") 
cc1.p <- cor(a[,i],a[,'GL'],method="pearson") 
empcol1.s <- paste(i,'.g.emp.spearman',sep='')
emp_rs[empRow,empcol1.s] <- cc1.s
empcol1.p <- paste(i,'.g.emp.pearson',sep='')
emp_rs[empRow,empcol1.p] <- cc1.p
i <- "GL"
empcol2 <- paste(i,'.b.emp.spearman',sep='')
emp_rs[empRow,empcol2] <- cc1.s
empcol1.p <- paste(i,'.b.emp.pearson',sep='')
emp_rs[empRow,empcol1.p] <- cc1.p

# Add rows for quantile cut offs 

emp_rs$replicate <- as.character(emp_rs$replicate)
emp_rs[2,1] <- "p0.025"
emp_rs[3,1] <- "p0.975"

# add gl and bf 

emp_rs[2,'GL.b.emp.spearman'] <- quantile(null_rs$BF.g.null.spearman, c(0.025, 0.975))[[1]]
emp_rs[3,'GL.b.emp.spearman'] <- quantile(null_rs$BF.g.null.spearman, c(0.025, 0.975))[[2]]
emp_rs[2,'BF.g.emp.spearman'] <- quantile(null_rs$GL.b.null.spearman, c(0.025, 0.975))[[1]]
emp_rs[3,'BF.g.emp.spearman'] <- quantile(null_rs$GL.b.null.spearman, c(0.025, 0.975))[[2]]

emp_rs[2,'GL.b.emp.pearson'] <- quantile(null_rs$BF.g.null.pearson, c(0.025, 0.975))[[1]]
emp_rs[3,'GL.b.emp.pearson'] <- quantile(null_rs$BF.g.null.pearson, c(0.025, 0.975))[[2]]
emp_rs[2,'BF.g.emp.pearson'] <- quantile(null_rs$GL.b.null.pearson, c(0.025, 0.975))[[1]]
emp_rs[3,'BF.g.emp.pearson'] <- quantile(null_rs$GL.b.null.pearson, c(0.025, 0.975))[[2]]


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


# Save each hypothesis to output in xlxs sheet for each dataset 

assign(h, sig2tp05)


# #ctrl shift c - to uncomment. then run at end
# #https://stackoverflow.com/questions/48799959/r-reduce-with-merge-and-more-than-2-suffixes-or-how-to-merge-multiple-datafr
# #Merge all datasets
# 
# list.df <- list(TvS,AIvSA,AIvSI,SAvSI)
# sfx <- c("_TvS", "_AIvSA", "_AIvSI", "_SAvSI")
# res <- list.df[[1]]
# for(i in head(seq_along(list.df), -1)) {
# 
#   res <- merge(res, list.df[[i+1]], all = TRUE,
#                suffixes = sfx[i:(i+1)], by = c("metaData"))
# }
# 
# # output excel sheet
# 
# library(xlsx)
# 
# fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_correlation.xlsx",sep='')
# 
# write.xlsx2(as.data.frame(res), file=fname, sheetName="ALL", row.names=FALSE)
# write.xlsx2(as.data.frame(TvS), file=fname, sheetName="TvS", append=TRUE, row.names=FALSE)
# write.xlsx2(as.data.frame(AIvSA), file=fname, sheetName="AIvSA", append=TRUE, row.names=FALSE)
# write.xlsx2(as.data.frame(AIvSI), file=fname, sheetName="AIvSI", append=TRUE, row.names=FALSE)
# write.xlsx2(as.data.frame(SAvSI), file=fname, sheetName="SAvSI", append=TRUE, row.names=FALSE)
# 
