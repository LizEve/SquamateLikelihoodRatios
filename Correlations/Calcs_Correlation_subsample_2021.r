
rm(list=ls())
# 2021 - do not have to redo this with fixed subsample
# 2021 - redoing with new burbrink and reeder 
##################################################################################################
## Test correlation across all hypotheses 
##################################################################################################

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh2021.RData")

# remove all but full dataset 
rm(list=setdiff(ls(), "all"))

##################################################################################################
## Test correlation across all hypotheses 
##################################################################################################
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# https://rpubs.com/aaronsc32/spearman-rank-correlation
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library(xlsx)
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
  if(length(a$Locus) < 5){
    c.d <- emptyFrame(h,s)
    return(c.d)
  }
  
  # Define which column numbers to test for correlation
  
  b <- match('BF',names(a))
  g <- match('GL',names(a))

  # Permute BF values in 1,000 separate columns, could permute BF or GL doesnt matter. 
  
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
## Subsample
##################################################################################################

all <- all %>% filter(dataSet == "Burbrink")

# full dataset numbers 

tvsx <- corTest(all,"TvS","all") 
aivsax <- corTest(all,"AIvSA","all")
aivsix <- corTest(all,"AIvSI","all")
savsix <- corTest(all,"SAvSI","all")
results <- bind_rows(tvsx,aivsax,aivsix,savsix) %>% filter(TestType == "emp.spearman")

r <- results

# run this x number of times 
reps = 1000

# subsample number and make sure variables are grouped
all <- all %>% group_by(variable)
sn = "150" # the name for the subsample - make sure to edit how many go into the subsample
print(sn)
for (n in 1:reps) {
  sb <- sample_n(all, 150)
  tvs <- corTest(sb,"TvS",sn) 
  aivsa <- corTest(sb,"AIvSA",sn)
  aivsi <- corTest(sb,"AIvSI",sn)
  savsi <- corTest(sb,"SAvSI",sn)
  results <- bind_rows(results,tvs,aivsa,aivsi,savsi) %>% filter(TestType == "emp.spearman")
  print(n)
}  

save(results,file="corsubsampleBurbrink")


all <- all %>% group_by(variable)
sn = "200"
print(sn)
for (n in 1:reps) {
  sb <- sample_n(all, 200)
  tvs <- corTest(sb,"TvS",sn) 
  aivsa <- corTest(sb,"AIvSA",sn)
  aivsi <- corTest(sb,"AIvSI",sn)
  savsi <- corTest(sb,"SAvSI",sn)
  results <- bind_rows(results,tvs,aivsa,aivsi,savsi) %>% filter(TestType == "emp.spearman")
  print(n)
}  

save(results,file="corsubsampleBurbrink")


all <- all %>% group_by(variable)
sn = "100"
print(sn)
for (n in 1:reps) {
  sb <- sample_n(all, 100)
  tvs <- corTest(sb,"TvS",sn) 
  aivsa <- corTest(sb,"AIvSA",sn)
  aivsi <- corTest(sb,"AIvSI",sn)
  savsi <- corTest(sb,"SAvSI",sn)
  results <- bind_rows(results,tvs,aivsa,aivsi,savsi) %>% filter(TestType == "emp.spearman")
  print(n)
}  

save(results,file="corsubsampleBurbrink")


##################################################################################################

all <- all %>% filter(dataSet == "Reeder")

# full dataset numbers 

tvsx <- corTest(all,"TvS","all") 
aivsax <- corTest(all,"AIvSA","all")
aivsix <- corTest(all,"AIvSI","all")
savsix <- corTest(all,"SAvSI","all")
results <- bind_rows(tvsx,aivsax,aivsix,savsix) %>% filter(TestType == "emp.spearman")

r <- results

# run this x number of times 
reps = 1000

# subsample number and make sure variables are grouped
all <- all %>% group_by(variable)
sn = "30" # the name for the subsample - make sure to edit how many go into the subsample
print(sn)
for (n in 1:reps) {
  sb <- sample_n(all, 30)
  tvs <- corTest(sb,"TvS",sn) 
  aivsa <- corTest(sb,"AIvSA",sn)
  aivsi <- corTest(sb,"AIvSI",sn)
  savsi <- corTest(sb,"SAvSI",sn)
  results <- bind_rows(results,tvs,aivsa,aivsi,savsi) %>% filter(TestType == "emp.spearman")
  print(n)
}  

save(results,file="corsubsampleReeder")


all <- all %>% group_by(variable)
sn = "25"
print(sn)
for (n in 1:reps) {
  sb <- sample_n(all, 25)
  tvs <- corTest(sb,"TvS",sn) 
  aivsa <- corTest(sb,"AIvSA",sn)
  aivsi <- corTest(sb,"AIvSI",sn)
  savsi <- corTest(sb,"SAvSI",sn)
  results <- bind_rows(results,tvs,aivsa,aivsi,savsi) %>% filter(TestType == "emp.spearman")
  print(n)
}  

save(results,file="corsubsampleReeder")


all <- all %>% group_by(variable)
sn = "35"
print(sn)
for (n in 1:reps) {
  sb <- sample_n(all, 35)
  tvs <- corTest(sb,"TvS",sn) 
  aivsa <- corTest(sb,"AIvSA",sn)
  aivsi <- corTest(sb,"AIvSI",sn)
  savsi <- corTest(sb,"SAvSI",sn)
  results <- bind_rows(results,tvs,aivsa,aivsi,savsi) %>% filter(TestType == "emp.spearman")
  print(n)
}  

save(results,file="corsubsampleReeder")



##################################################################################################
## Graph
##################################################################################################
rm(list=ls())

mytheme <- theme_bw() + theme(panel.border = element_blank()) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(), # get rid of major grid
    plot.title = element_text(hjust = 0.5),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.background = element_rect(colour = "transparent",fill = "transparent"))

library("readxl")
library(dplyr)
library(ggplot2)
library("viridis")

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("corsubsample")
setwd("/Users/chatnoir/Projects/Squam/Graphs/")
empvals <- read_excel("/Users/chatnoir/Projects/Squam/Graphs/Tables/AllSchmoosh.xlsx")

all <- bind_rows(results,empvals) %>% filter(Sample != "50l") %>% filter(Sample != "50h") %>%
  mutate(Sample=factor(Sample, levels=c("100","500","1k","2k","3k","4k","sig50m","sig60m","sig80m","sig0x","sig1x","sig2x","all")))

tvs <- all %>% filter(Hypothesis == "TvS")
  
df <- all %>% filter(Hypothesis == "AIvSA")

#df <- all

x.tic <- seq(-lowerx,upperx,ticx) 
y.tic <- seq(-lowery,uppery,ticy)
#x.tic <- c(-lowerx,upperx,ticx) 
#y.tic <- c(-lowery,uppery,ticy)


scat <- ggplot(df, aes(x=Loci,y=corr.eff)) + 
  geom_point(alpha=0.75,aes(color=as.factor(Sample))) + mytheme

scat
  
  
scat <- ggplot(df, aes(x=Loci,y=corr.eff)) + 
  geom_point(alpha=0.5,aes(color=p)) +
  scale_color_viridis(direction = -1,option='D') + mytheme

scat
  
   +
  scale_color_viridis(direction = -1,option='D')
  #coord_cartesian(ylim=y.tic,xlim = x.tic) +
  #scale_y_continuous(breaks = y.tic) + 
  #scale_x_continuous(breaks = x.tic) +
  #labs(x="dGLS",y="ln(BF)", color="Metric Disagreement") + 
  #scale_color_manual(values=cc) 
#quartz()
scat


