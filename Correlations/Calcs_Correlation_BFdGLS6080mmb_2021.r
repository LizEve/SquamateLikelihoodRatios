rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library(xlsx)

##################################################################################################
## Test correlation across all hypotheses 
##################################################################################################
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# https://rpubs.com/aaronsc32/spearman-rank-correlation

emptyFrame <- function(h,s){
  c.d <- data.frame("test"=c("emp.spearman","emp.pearson"),
                    "sample"=c(s,s),
                    "hypothesis"=c(h,h),
                    "Correlation"=c(NA,NA),
                    "emp"=c(NA,NA),
                    "p0.975"=c(NA,NA),
                    "p0.95"=c(NA,NA),
                    "p0.05"=c(NA,NA),
                    "p0.025"=c(NA,NA))
  return(c.d)
}

corTest <- function(d,h,s){
  
  # Create two datasets, one to add permutation columns to, one to pull from. because that's how I coded it?
  a <- d %>% filter(variable == h) 
  perms <- d %>% filter(variable == h) 
  nloci <- length(a$Locus)
  maxB <- max(a$BF)
  minB <- min(a$BF)
  maxG <- max(a$GL)
  minG <- min(a$GL)
  #print(length(a$Locus))
  
  # Safety check for if a subset of data has too few observations
  if(length(a$Locus) < 1){
    c.d <- emptyFrame(h,s)
    return(c.d)
  }
  
  # Define which column numbers to test for correlation
  
  b <- match('BF',names(a))
  g <- match('GL',names(a))
  
  # Permute BF values in 1,000 seperate columns, could permute BF or GL doesnt matter. 
  
  for (n in 1:1000) {
    newcolname1 <- paste('BF',n,sep='') 
    perms[, newcolname1] <- sample(a[[b]], replace = FALSE)
    #newcolname2 <- paste('GL',n,sep='') 
    #null[, newcolname2] <- sample(a[[g]], replace = FALSE)
  }
  
  # Create new dataset to add null rank coefficients to 
  
  assign('null', setNames(data.frame(1:1000), "replicate"))
  
  # Add bf gl 
  
  for (n in 1:1000) {
    
    # Grab column names, calculate value 
    colname1 <- paste('BF',n,sep='') 
    cc1.s <- cor(perms[,'GL'],perms[,colname1],method="spearman") 
    null[n,'null.spearman'] <- cc1.s
    cc1.p <- cor(perms[,'GL'],perms[,colname1],method="pearson") 
    null[n,'null.pearson'] <- cc1.p
  }
  
  
  # Dataframe for emp scores
  emp <- data.frame(TestType=character(),
                    corr.eff=integer(),
                    p=integer(),
                    p0.975=integer(),
                    p0.95=integer(),
                    p0.05=integer(),
                    p0.025=integer(),stringsAsFactors=FALSE)
  
  # Designate which row is for which test 
  
  emp[1,'TestType'] <- 'emp.spearman'
  emp[2,'TestType'] <- 'emp.pearson'
  
  
  # Calculate correlation stats from emperical data,  add to df
  
  e.s <- cor(a[[b]],a[[g]],method="spearman") 
  e.p <- cor(a[[b]],a[[g]],method="pearson") 
  
  emp[1,'corr.eff'] <- e.s
  emp[2,'corr.eff'] <- e.p
  
  
  # Calculate 2 sided p values
  
  p.s <- sum(abs(null$null.spearman) >= abs(e.s)) / length(null$null.spearman)
  p.p <- sum(abs(null$null.pearson) >= abs(e.s)) / length(null$null.pearson)
  
  emp[1,'p'] <- p.s
  emp[2,'p'] <- p.p
  
  
  # Add rows for quantile cut offs 
  
  emp[1,'p0.025'] <- quantile(null$null.spearman, c(0.025, 0.975),na.rm = TRUE)[[1]]
  emp[1,'p0.975'] <- quantile(null$null.spearman, c(0.025, 0.975),na.rm = TRUE)[[2]]
  emp[1,'p0.05'] <- quantile(null$null.spearman, c(0.05, 0.95),na.rm = TRUE)[[1]]
  emp[1,'p0.95'] <- quantile(null$null.spearman, c(0.05, 0.95),na.rm = TRUE)[[2]]
  emp[2,'p0.025'] <- quantile(null$null.pearson, c(0.025, 0.975),na.rm = TRUE)[[1]]
  emp[2,'p0.975'] <- quantile(null$null.pearson, c(0.025, 0.975),na.rm = TRUE)[[2]]
  emp[2,'p0.05'] <- quantile(null$null.pearson, c(0.05, 0.95),na.rm = TRUE)[[1]]
  emp[2,'p0.95'] <- quantile(null$null.pearson, c(0.05, 0.95),na.rm = TRUE)[[2]]
  
  # Add significance yes/no and other data 
  
  corr.data <- emp %>%
    mutate_if(is.numeric,round,2) %>% 
    add_column(Hypothesis = h, .after = "TestType") %>% 
    add_column(Sample = s, .after = "TestType") %>% 
    add_column(Loci = nloci, .after = "p") %>%  
    add_column(minG = minG, .after = "Loci")%>% 
    add_column(maxG = maxG, .after = "Loci")%>% 
    add_column(minB = minB, .after = "Loci")%>% 
    add_column(maxB = maxB, .after = "Loci")
  
  return(corr.data)
}

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


# Originally ran with threshold of dGLS being 0.5, 2x = 1, 4x=2
# trying again with 2, 2x=4, 4x=8

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
## run tests 
##################################################################################################


tvs <- corTest(all.exp,"TvS","all") 
aivsa <- corTest(all.exp,"AIvSA","all")
aivsi <- corTest(all.exp,"AIvSI","all")
savsi <- corTest(all.exp,"SAvSI","all")




tvs.60m <- corTest(all.60m,"TvS","sig60m")
aivsa.60m <- corTest(all.60m,"AIvSA","sig60m")
aivsi.60m <- corTest(all.60m,"AIvSI","sig60m")
savsi.60m <- corTest(all.60m,"SAvSI","sig60m")


tvs.80m <- corTest(all.80m,"TvS","sig80m")
aivsa.80m <- corTest(all.80m,"AIvSA","sig80m")
aivsi.80m <- corTest(all.80m,"AIvSI","sig80m")
savsi.80m <- corTest(all.80m,"SAvSI","sig80m")


tvs.60mb <- corTest(all.60mb,"TvS","sig60mb")
aivsa.60mb <- corTest(all.60mb,"AIvSA","sig60mb") # streicher only has 0 for GL for 60mb, so it shows up as NA
aivsi.60mb <- corTest(all.60mb,"AIvSI","sig60mb")
savsi.60mb <- corTest(all.60mb,"SAvSI","sig60mb")




tvs.80mb <- corTest(all.80mb,"TvS","sig80mb")
aivsa.80mb <- corTest(all.80mb,"AIvSA","sig80mb")
aivsi.80mb <- corTest(all.80mb,"AIvSI","sig80mb")
savsi.80mb <- corTest(all.80mb,"SAvSI","sig80mb")




results <- bind_rows(tvs,aivsa,aivsi,savsi,
                        tvs.80m,aivsa.80m,aivsi.80m,savsi.80m,
                        tvs.60m,aivsa.60m,aivsi.60m,savsi.60m,
                        tvs.80mb,aivsa.80mb,aivsi.80mb,savsi.80mb,
                        tvs.60mb,aivsa.60mb,aivsi.60mb,savsi.60mb) %>% filter(TestType == "emp.spearman")


# output excel sheet

library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/Tables/Correlation_",dataset,"_mb_2021.xlsx",sep='')

write.xlsx2(as.data.frame(results), file=fname, sheetName=dataset, row.names=FALSE,append=TRUE)

