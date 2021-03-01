rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library("viridis")
# set dataset 

#dataset <- "SinghalOG"
dataset <- "Streicher"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Colors 
# https://www.w3schools.com/colors/colors_picker.asp
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS"))
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))

# Get total number of loci 
loci <- length(mLGL$Locus)

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10 - not including average calcs
names(mLBF[7:10])
mLBF[7:10] <- mLBF[7:10]/2
# for singhal OG 
#names(mLBF[4])
#mLBF[4] <- mLBF[4]/2


# Want to stack datasets so GL and BF columns have both comparisons for each hypothesis. 

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
AI.g <- bind_rows(reformatdf(mLGL,"AI","AIvSA",1,c(1,13)),
                  reformatdf(mLGL,"AI","AIvSI",1,c(1,13))) %>% select(-SIvAvg)

SA.g <- bind_rows(reformatdf(mLGL,"SA","AIvSA",-1,c(1,13)),
                  reformatdf(mLGL,"SA","SAvSI",1,c(1,13))) %>% select(-SIvAvg)

SI.g <- bind_rows(reformatdf(mLGL,"SI","AIvSI",-1,c(1,13)),
                  reformatdf(mLGL,"SI","SAvSI",-1,c(1,13))) %>% select(-SIvAvg)

AI.b <- bind_rows(reformatdf(mLBF,"AI","AIvSA",1,c(1,13)),
                  reformatdf(mLBF,"AI","AIvSI",1,c(1,13))) %>% select(-SIvAvg)

SA.b <- bind_rows(reformatdf(mLBF,"SA","AIvSA",-1,c(1,13)),
                  reformatdf(mLBF,"SA","SAvSI",1,c(1,13))) %>% select(-SIvAvg)

SI.b <- bind_rows(reformatdf(mLBF,"SI","AIvSI",-1,c(1,13)),
                  reformatdf(mLBF,"SI","SAvSI",-1,c(1,13))) %>% select(-SIvAvg)

TS.g <- reformatdf(mLGL,"TS","TvS",1,c(1,13)) %>% select(-SIvAvg)
TS.b <- reformatdf(mLBF,"TS","TvS",1,c(1,13)) %>% select(-SIvAvg)

# Singhal OG
#TS.g <- reformatdf(mLGL,"TS","TvS",1,c(1,3))
#TS.b <- reformatdf(mLBF,"TS","TvS",1,c(1,3))

# MERGE 

AI.a <- left_join(AI.g, AI.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
SA.a <- left_join(SA.g, SA.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
SI.a <- left_join(SI.g, SI.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
TS.a <- left_join(TS.g, TS.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))

rm(AI.g,AI.b,SA.b,SA.g,SI.b,SI.g,TS.b,TS.g)

#----COMBINE-----------------------------------------------------------------------------------------------------------------------------------------------------

# Create dataframe with all metadata, locus names, one col for hypothesis, one call for bf, and one gl
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

zeroAI <- AI.a %>% filter(AI.g==0)
zeroSA <- SA.a %>% filter(SA.g==0)
zeroSI <- SI.a %>% filter(SI.g==0)
#zeroG <- AI.a %>% filter(between(AI.g,-0.01,0.01))

#----Trees-----------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(ape)
library(phangorn)
library(phytools)
library(geiger)

makeTree <- function(c.tree,title,keep){
  y <- ladderize(c.tree, right = TRUE)
  y.root <- root(y,"gallus_gallus")
  pr.tree<-drop.tip(y.root,
                    setdiff(y.root$tip.label,keep))
  plot(pr.tree,cex=0.4,no.margin=TRUE,
       show.node.label = TRUE,
       use.edge.length = FALSE,
       node.depth = 2,
       edge.width = 0.5)
  edgelabels(round(pr.tree$edge.length,7), bg="white", 
             col="black", font=0.4, cex=0.5,  
             adj = c(0.5, -0.5))
  title(main=title,cex.main = 0.5, font.main = 3, line=-2)
}

makePdf <- function(fname,c){
  fpath <- paste('/Users/ChatNoir/Projects/Squam/Streicher/BranchLenInvestigation/',fname,sep = '')
  f <- readLines(fpath)
  g <- grep('Tree in newick format:',f)+2
  t <- read.tree(text=f[g])
  pdf(file=paste(c,".pdf", sep=""),width = 3, height = 4)
  makeTree(t,c, k)
  dev.off()
}

dataset <- "Streicher"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))


k=c('gonatodes_sp','anniella_pulchra', 'heloderma_suspectum', 'varanus_exanthematicus', 'xenosaurus_platyceps', 'lanthanotus_borneensis', 'hydrosaurus_sp', 'anolis_carolinensis_2', 'uta_stansburiana', 'micrurus_fulvius', 'python_molurus_2', 'typhlops_jamaicensis')

#quartz()
f <- 'uce-1758.toxicoferaAI.constraint.iqtree'
makePdf(f,"IQ_1758_AI")
f <- 'uce-1758.toxicoferaSA.constraint.iqtree'
makePdf(f,"IQ_1758_SA")

f <- 'uce-7785.toxicoferaAI.constraint.iqtree'
makePdf(f,"IQ_7785_AI")
f <- 'uce-7785.toxicoferaSI.constraint.iqtree'
makePdf(f,"IQ_7785_SI")

#----SS-Trees-----------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(ape)
library(phangorn)
library(phytools)
library(geiger)
dataset <- "Streicher"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

k=c('gonatodes_sp','anniella_pulchra', 'heloderma_suspectum', 'varanus_exanthematicus', 'xenosaurus_platyceps', 'lanthanotus_borneensis', 'hydrosaurus_sp', 'anolis_carolinensis_2', 'uta_stansburiana', 'micrurus_fulvius', 'python_molurus_2', 'typhlops_jamaicensis')
#quartz()

# take last 70 trees, they should be closest to posterior. Log file said it went from prior to posterior, but first tree had equal branche lenghts, so I think the trees are written to file from prior to posterior
makeTree <- function(c.tree,title,keep){
  y <- ladderize(c.tree, right = TRUE)
  y.root <- root(y,"gallus_gallus")
  pr.tree<-drop.tip(y.root,
                    setdiff(y.root$tip.label,keep))
  plot(pr.tree,cex=0.4,no.margin=TRUE,
       show.node.label = TRUE,
       use.edge.length = FALSE,
       node.depth = 2,
       edge.width = 0.5)
  edgelabels(round(pr.tree$edge.length,4), bg="white", 
             col="black", font=0.4, cex=0.5,  
             adj = c(0.5, -0.5))
  title(main=title,cex.main = 0.5, font.main = 3, line=-2)
}

pdfMake <- function(fname,cap,keep){
  fpath <- paste('/Users/ChatNoir/Projects/Squam/Streicher/BranchLenInvestigation/StreicherSS/',fname,sep = '')
  ts <- read.nexus(fpath)
  trees <- ts[3900:3975]
  c.t <- consensus.edges(trees,p=0.5)
  pdf(file=paste(cap,"_prior.pdf", sep=""),width = 7, height = 3)
  par(mfrow=(c(1,5)))
  makeTree(c.t,"prior consensus",keep)
  makeTree(ts[[3975]],paste("3975, ",cap,sep=""),keep)
  makeTree(ts[[3950]],paste("3950, ",cap,sep=""),keep)
  makeTree(ts[[3920]],paste("3920, ",cap,sep=""),keep)
  makeTree(ts[[3900]],paste("3900, ",cap,sep=""),keep)
  dev.off()
  # Get consensus for 70 near prior
  # Get tree 70,60,50,40,30,20,10
  trees <- ts[2:75]
  c.t <- consensus.edges(trees,p=0.5)
  pdf(file=paste(cap,"_posterior.pdf", sep=""),width = 7, height = 3)
  par(mfrow=(c(1,5)))
  makeTree(c.t,"posterior consensus",keep)
  makeTree(ts[[2]],paste("0, ",cap,sep=""),keep)
  makeTree(ts[[20]],paste("20, ",cap,sep=""),keep)
  makeTree(ts[[50]],paste("50, ",cap,sep=""),keep)
  makeTree(ts[[70]],paste("70, ",cap,sep=""),keep)
  dev.off()
}


f <- 'uce-1758/uce-1758.toxai.run1.t'
pdfMake(f,"AI_1758",k) 
f <- 'uce-1758/uce-1758.toxsa.run1.t'
pdfMake(f,"SA_1758",k) 

f <- 'uce-7785/uce-7785.toxai.run1.t'
pdfMake(f,"AI_7785",k) 
f <- 'uce-7785/uce-7785.toxsi.run1.t'
pdfMake(f,"SI_7785",k) 



#----SS-Trees-Working Notes----------------------------------------------------------------------------------------------------------------------------------------------------
y <- ladderize(t, right = TRUE)
y.root <- root(y,"gallus_gallus")
pr.tree<-drop.tip(y.root,
                  setdiff(y.root$tip.label,keep))

plot(pr.tree,cex=0.5,no.margin=TRUE,show.node.label = TRUE)

edgelabels(round(pr.tree$edge.length,8), bg="white", 
           col="black", font=0.5, cex=0.5, frame='none', 
           adj = c(0.5, 0.0))



trees <- ts[3900:3975]
trees <- ts[0:70]
c.tree <- consensus.edges(trees,p=0.5)
#quartz()
#c.tree <- consensus.edges(trees,p=0.3)
y <- ladderize(c.tree, right = TRUE)
y <- ladderize(ts[[50]], right = TRUE)
quartz()
y.root <- root(y,"gallus_gallus")
pr.tree<-drop.tip(y.root,
                  setdiff(y.root$tip.label,keep))
pr.tree <- y.root

plot(pr.tree,cex=0.4,no.margin=FALSE,
     show.node.label = TRUE,node.depth = 2,
     use.edge.length = FALSE,
     edge.width = 0.5)

edgelabels(round(pr.tree$edge.length,5), bg="white", 
           col="black", font=0.5, cex=0.5,  
           adj = c(0.5, -0.5))
title(main="title",cex.main = 0.5, font.main = 3)



# Get consensus for 70 near prior
# Get tree 70,50,20,0
#quartz()
trees <- ts[3900:3975]
c.t <- consensus.edges(trees,p=0.5)
pdf(file="AI_1758_prior.pdf",width = 7, height = 3)
par(mfrow=(c(1,5)))
makeTree(c.t,"prior consensus")
makeTree(ts[[3975]],"3975, AI 1758")
makeTree(ts[[3950]],"3950, AI 1758")
makeTree(ts[[3920]],"3920, AI 1758")
makeTree(ts[[3900]],"3900, AI 1758")
dev.off()
# Get consensus for 70 near prior
# Get tree 70,60,50,40,30,20,10
trees <- ts[0:70]
c.t <- consensus.edges(trees,p=0.5)
pdf(file="AI_1758_posterior.pdf",width = 7, height = 3)
par(mfrow=(c(1,5)))
makeTree(c.t,"posterior consensus")
makeTree(ts[[2]],"0, AI 1758")
makeTree(ts[[20]],"20, AI 1758")
makeTree(ts[[50]],"50, AI 1758")
makeTree(ts[[70]],"70, AI 1758")
dev.off()