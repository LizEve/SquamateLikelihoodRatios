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
  maxB <- max(a$drank.b)
  minB <- min(a$drank.b)
  maxG <- max(a$drank.g)
  minG <- min(a$drank.g)
  #print(length(a$Locus))
  
  # Safety check for if a subset of data has too few observations
  if(length(a$Locus) < 5){
    c.d <- emptyFrame(h,s)
    return(c.d)
  }
  
  # Define which column numbers to test for correlation
  
  b <- match('drank.b',names(a))
  g <- match('drank.g',names(a))
  
  # Permute BF values in 1,000 separate columns, could permute BF or GL doesnt matter. 
  
  for (n in 1:1000) {
    newcolname1 <- paste('drank.b',n,sep='') 
    perms[, newcolname1] <- sample(a[[b]], replace = FALSE)
    #newcolname2 <- paste('GL',n,sep='') 
    #null[, newcolname2] <- sample(a[[g]], replace = FALSE)
  }
  
  # Create new dataset to add null rank coefficients to 
  
  assign('null', setNames(data.frame(1:1000), "replicate"))
  
  # Add bf gl 
  
  for (n in 1:1000) {
    
    # Grab column names, calculate value 
    colname1 <- paste('drank.b',n,sep='') 
    cc1.s <- cor(perms[,'drank.g'],perms[,colname1],method="spearman") 
    null[n,'null.spearman'] <- cc1.s
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
  
  
  
  # Calculate correlation stats from emperical data,  add to df
  
  e.s <- cor(a[[b]],a[[g]],method="spearman") 
  
  emp[1,'corr.eff'] <- e.s
  
  # Calculate 2 sided p values
  
  p.s <- sum(abs(null$null.spearman) >= abs(e.s)) / length(null$null.spearman)
  
  emp[1,'p'] <- p.s
  
  # Add rows for quantile cut offs 
  
  emp[1,'p0.025'] <- quantile(null$null.spearman, c(0.025, 0.975),na.rm = TRUE)[[1]]
  emp[1,'p0.975'] <- quantile(null$null.spearman, c(0.025, 0.975),na.rm = TRUE)[[2]]
  emp[1,'p0.05'] <- quantile(null$null.spearman, c(0.05, 0.95),na.rm = TRUE)[[1]]
  emp[1,'p0.95'] <- quantile(null$null.spearman, c(0.05, 0.95),na.rm = TRUE)[[2]]
  
  
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
## get dataset 
##################################################################################################

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh2021April.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/")


all <- all %>% filter(dataSet == "Streicher")

# full dataset numbers 

tvsx <- corTest(all,"TvS","all") 
aivsax <- corTest(all,"AIvSA","all")
aivsix <- corTest(all,"AIvSI","all")
savsix <- corTest(all,"SAvSI","all")
results <- bind_rows(tvsx,aivsax,aivsix,savsix) %>% filter(TestType == "emp.spearman")

r <- results


##################################################################################################
## linear regression
##################################################################################################


rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library(xlsx)
setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh2021April.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/")

d <- all %>% filter(dataSet == "Burbrink")


#y-gl
#x-bf
lr.df<- data.frame(hypothesis=character(),
                   subset=character(),
                   a=numeric(), 
                   b=numeric(), 
                   r2=numeric(),
                   nloci=numeric(),
                   source=character())
hypos <- c("TvS","AIvSA","AIvSI","SAvSI")


v <- "variable"
names(d)



for (h in hypos){
  print(h)
  df <- d %>% filter(variable == h)
  b <- "drank.b"
  print(max(df[[b]]))
  g <- "drank.g"
  m = lm(df[[b]] ~ df[[g]])
  print(m)
  print(summary(m)$r.squared)
  l <- list(hypothesis=h, 
            a = format(unname(coef(m)[1]), digits = 2),
            b = format(unname(coef(m)[2]), digits = 2), 
            r2 = format(summary(m)$r.squared, digits = 3),
            nloci = length(df[[g]]),
            source = dataset)
  #print(data.frame(l))
  lr.df <- rbind(lr.df,data.frame(l))
}



fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/Tables/RankLinearRegression_Burbrink_2021.xlsx",sep='')

write.xlsx2(as.data.frame(lr.df), file=fname, sheetName=dataset, row.names=FALSE,append=TRUE)
