rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

writeCSV <- function (outData,outFile){
  write.csv(outData, file=outFile, row.names=FALSE)
}

BFofAvg <- function(numeratorConstraint,margLikes){
  
  # Rename column name for numberator in BF calculation. 
  colnames(margLikes)[colnames(margLikes)==numeratorConstraint] <- "numeratorConstraint"
  
  # Start a list of BF values for all loci in dataset
  numLoci <- length(margLikes$numeratorConstraint)
  
  # Make lists and add to data frame with locus names, constraint name, and 0 BFs
  numeratorList <- rep(c(numeratorConstraint), times=numLoci)
  outBFList <- vector(mode="numeric",length=numLoci)
  
  # Make dataframe 
  outBFdf <- subset(margLikes, select = c('Locus'))
  outBFdf$Hypothesis <- numeratorList
  outBFdf$BF <- outBFList
    
  # create a dataset of alternate constraints only to average over
  drops <- c("Locus","numeratorConstraint")
  df <- margLikes[,!(names(margLikes) %in% drops)]
  
  # for locus in dataset 
  for (i in 1:nrow(outBFdf)){
    
    # Get all alternate marginal likehoods 
    otherMargLikes <- as.numeric(df[i,])
    
    # Maximum of alternate margLH
    a <- max(otherMargLikes)
    
    # Calculate average across alternate constraints
    avgOtherMargLikes <- a+log(sum(exp(otherMargLikes-a)))-log(length(names(df)))
    
    # Calculate bayes factor - 2ln(BF)
    outBFdf[i,3] <- 2*(margLikes$numeratorConstraint[i] - avgOtherMargLikes)}
  
  return(outBFdf)
}



setwd("/Users/ChatNoir/Projects/Squam/Graphs/BF")
mLBF <- read.table("StreicherMargLikeSummary.txt",header=TRUE)
mLdGLS <- read.table("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/StreicherdGLSSummary.txt",header=TRUE)
length(mLBF$Locus[!(mLBF$Locus %in% mLdGLS$Locus)])
#same lenght all good. 
mL <- read.table("StreicherMargLikeSummary.txt",header=TRUE)


#mLBF <- read.table("SquamNTGMargLikeSummary.txt",header=TRUE)
#mLdGLS <- read.table("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/SquamNTGdGLSSummary.txt",header=TRUE)
#length(mLBF$Locus[!(mLBF$Locus %in% mLdGLS$Locus)])
#removeLoci <- mLBF$Locus[!(mLBF$Locus %in% mLdGLS$Locus)]
#mLa <- read.table("SquamNTGMargLikeSummary.txt",header=TRUE)
#mL <- mLa[! mLa$Locus %in% removeLoci,]

# Read in table of marginal Likelihoods

mLTox <- mL[,!(names(mL) %in% c("Sclero"))]
names(mL)
# dataframes with only sclero and one tox hypoth 
mLaiS <- mL[,(names(mL) %in% c("ToxAngIg","Sclero","Locus"))]
mLsaS <- mL[,(names(mL) %in% c("ToxSnAng","Sclero","Locus"))]
mLsiS <- mL[,(names(mL) %in% c("ToxSnIg","Sclero","Locus"))]



# Get list of BF for each relationship, using an average of all other constraints as a composite hypothesis. 

# Pairwise with sclero 
saiBF <- subset(mLaiS, select = c('Locus'))
saiBF$Hypothesis <- rep(c("SvAI"), times=length(mLaiS$Locus))
saiBF$BF <- 2*(mLaiS$Sclero-mLaiS$ToxAngIg) 

ssaBF <- subset(mLsaS, select = c('Locus'))
ssaBF$Hypothesis <- rep(c("SvSA"), times=length(mLsaS$Locus))
ssaBF$BF <- 2*(mLsaS$Sclero-mLsaS$ToxSnAng) 

ssiBF <- subset(mLsiS, select = c('Locus'))
ssiBF$Hypothesis <- rep(c("SvSI"), times=length(mLsiS$Locus))
ssiBF$BF <- 2*(mLsiS$Sclero-mLsiS$ToxSnIg) 

# all
ScleroBF <- BFofAvg("Sclero",mL)
ToxAngIgBF  <- BFofAvg("ToxAngIg",mL)
ToxSnAngBF  <- BFofAvg("ToxSnAng",mL) 
ToxSnIgBF <- BFofAvg("ToxSnIg",mL)

# differences only within tox constraints 
ToxAngIgBFtx  <- BFofAvg("ToxAngIg",mLTox)
ToxSnAngBFtx  <- BFofAvg("ToxSnAng",mLTox) 
ToxSnIgBFtx <- BFofAvg("ToxSnIg",mLTox)

# Combine and factor
allBF <- bind_rows(ScleroBF,ToxAngIgBF,ToxSnAngBF,ToxSnIgBF)
toxBF <- bind_rows(ToxAngIgBFtx,ToxSnAngBFtx,ToxSnIgBFtx)
allBF$Hypothesis <- factor(allBF$Hypothesis, levels=c("Sclero","ToxAngIg","ToxSnAng","ToxSnIg"))
toxBF$Hypothesis <- factor(toxBF$Hypothesis, levels=c("ToxAngIg","ToxSnAng","ToxSnIg"))
allBF <- allBF %>% separate(Locus, c("LocusType","LocusNumber"), "-", remove = FALSE)
toxBF <- toxBF %>% separate(Locus, c("LocusType","LocusNumber"), "-", remove = FALSE)

pairBF <- bind_rows(saiBF,ssaBF,ssiBF)

# Support counts 

allBF <- mutate(allBF, Support = case_when(BF >= 10 ~ "Strong",
                                           BF > -10 & BF < 10 ~ "Ambiguous",
                                           BF <= -10 ~ "Strong_Against")) %>%
  mutate(Percentile = ntile(desc(BF),100)) %>% mutate(locusRank=dense_rank(desc(BF)))
toxBF <- mutate(toxBF, Support = case_when(BF >= 10 ~ "Strong",
                                           BF > -10 & BF < 10 ~ "Ambiguous",
                                           BF <= -10 ~ "Strong_Against")) %>%
  mutate(Percentile = ntile(desc(BF),100)) %>% mutate(locusRank=dense_rank(desc(BF)))
pairBF <- mutate(pairBF, Support = case_when(BF >= 10 ~ "Strong",
                                           BF > -10 & BF < 10 ~ "Ambiguous",
                                           BF <= -10 ~ "Strong_Against")) %>%
  mutate(Percentile = ntile(desc(BF),100)) %>% mutate(locusRank=dense_rank(desc(BF)))
allSupport <- allBF %>% group_by(Hypothesis,Support) %>% dplyr::summarise(counts=n())
allSupport$Support <- factor(allSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))
toxSupport <- toxBF %>% group_by(Hypothesis,Support) %>% dplyr::summarise(counts=n())
toxSupport$Support <- factor(toxSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))
pairSupport <- pairBF %>% group_by(Hypothesis,Support) %>% dplyr::summarise(counts=n())
pairSupport$Support <- factor(pairSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))

# Save R data here 


# Absolute values, averages of differences 
margLHablsAve <- mL %>% mutate(absAvg=(abs(Sclero-ToxAngIg)+
                                         abs(Sclero-ToxSnAng)+
                                         abs(Sclero-ToxSnIg)+
                                         abs(ToxAngIg-ToxSnAng)+
                                         abs(ToxAngIg-ToxSnIg)+
                                         abs(ToxSnAng-ToxSnIg))/6) %>% select(Locus,absAvg)




quartz()


############## Saved Plots ############## ############## ############## ############## ############## 

# Streicher_BF_Violin
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
txList <- c("#38598CFF","#2BB07FFF","#C2DF23FF")
all <- ggplot(allBF, aes(x=Hypothesis, y=BF, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + scale_color_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank())

all
all.trunc <- ggplot(allBF, aes(x=Hypothesis, y=BF, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + ylim(-50,50) + stat_summary(fun.data=data_summary) + 
  scale_color_manual(values=myList) + theme_bw() + 
  theme_bw() + theme(panel.border = element_blank())
all.trunc

a <- ggarrange(all,all.trunc, ncol=2, nrow=1, common.legend = TRUE, legend="right") 
ggsave("Streicher_BF_Violin_all.jpg", plot=a, width = 40, height = 20, units = "cm")




tox <- ggplot(toxBF, aes(x=Hypothesis, y=BF, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + scale_color_manual(values=txList) + 
  theme_bw() + theme(panel.border = element_blank())
tox.trunc <- ggplot(toxBF, aes(x=Hypothesis, y=BF, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + ylim(-25,25) + stat_summary(fun.data=data_summary) + 
  scale_color_manual(values=txList) + 
  theme_bw() + theme(panel.border = element_blank())
tox
tox.trunc
t <- ggarrange(tox,tox.trunc, ncol=2, nrow=1, common.legend = TRUE, legend="right")

ggsave("Streicher_BF_Violin_tox.jpg", plot=t, width = 40, height = 20, units = "cm")

########### Streicher_BF_Histo ###########

myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
txList <- c("#38598CFF","#2BB07FFF","#C2DF23FF")

all <- ggplot(allBF, aes(x=BF, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=myList)+scale_fill_manual(values=myList)
all
all.trunc <- ggplot(allBF, aes(x=BF, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(0,10),color=c("grey55", "black"), linetype="dashed", size=0.5) + 
  xlim(-50,50) + scale_color_manual(values=myList)+scale_fill_manual(values=myList)

all.trunc 

a <- ggarrange(all,all.trunc, ncol=1, nrow=2, common.legend = TRUE, legend="right") 
a
ggsave("Streicher_BF_Histo_all.jpg", plot=a, width = 40, height = 20, units = "cm")


tox <- ggplot(toxBF, aes(x=BF, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=txList)+scale_fill_manual(values=txList) +
  theme_bw() + theme(panel.border = element_blank()) + 
  coord_cartesian(xlim=c(-90,90)) + scale_x_continuous(breaks=seq(-90,90,10))
tox.trunc <- ggplot(toxBF, aes(x=BF, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw() + theme(panel.border = element_blank()) + 
  geom_vline(xintercept=c(0,10),color=c("grey55", "black"), linetype="dashed", size=0.5) + 
  xlim(-25,25) + scale_color_manual(values=txList)+scale_fill_manual(values=txList)

t <- ggarrange(tox,tox.trunc, ncol=1, nrow=2, common.legend = TRUE, legend="right") 
t
ggsave("Streicher_BF_Histo_tox.jpg", plot=t, width = 40, height = 20, units = "cm")


# sclero avg vs pairwise 
sc <- ggplot(ScleroBF, aes(fill=Hypothesis, x=BF, color=Hypothesis)) + 
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(-10,10),color=c("black", "black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values="orange") + scale_fill_manual(values="orange") +
  theme_bw() + theme(panel.border = element_blank()) + 
  coord_cartesian(xlim=c(-110,100)) + scale_x_continuous(breaks=seq(-110,100,10)) +
  labs(x="BF Sclero vs avg")
sc
sai <- ggplot(SaiBF, aes(fill=Hypothesis, x=BF, color=Hypothesis)) + 
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(-10,10),color=c("black", "black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values="orange") + scale_fill_manual(values="orange") +
  theme_bw() + theme(panel.border = element_blank()) + 
  coord_cartesian(xlim=c(-110,100)) + scale_x_continuous(breaks=seq(-110,100,10)) +
  labs(x="BF Sclero vs AngIg")
ssa <- ggplot(SsaBF, aes(fill=Hypothesis, x=BF, color=Hypothesis)) + 
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(-10,10),color=c("black", "black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values="orange") + scale_fill_manual(values="orange") +
  theme_bw() + theme(panel.border = element_blank()) + 
  coord_cartesian(xlim=c(-110,100)) + scale_x_continuous(breaks=seq(-110,100,10)) +
  labs(x="BF Sclero vs SnAng")
ssi <- ggplot(SsiBF, aes(fill=Hypothesis, x=BF, color=Hypothesis)) + 
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(-10,10),color=c("black", "black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values="orange") + scale_fill_manual(values="orange") +
  theme_bw() + theme(panel.border = element_blank()) + 
  coord_cartesian(xlim=c(-110,100)) + scale_x_continuous(breaks=seq(-110,100,10)) +
  labs(x="BF Sclero vs SnIg")

a <- ggarrange(sc,sai,ssi,ssa, ncol=1, nrow=4, legend = "none")
a

#ggsave("Streicher_BF_Histo_ScleroPairwise.jpg", plot=a, width = 20, height = 20, units = "cm")




##########Streicher_BF_Classic###########

myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
txList <- c("#38598CFF","#2BB07FFF","#C2DF23FF")

SC <- ggplot(data=ScleroBF) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-BF,sum),y=BF, color = BF < 0, fill = BF < 0) +
  scale_color_manual(values=c("orange","grey40")) + scale_fill_manual(values=c("orange","grey40"))+ 
  labs(y="BF Hypoth - Average", title="Scleroglossa")
AI <- ggplot(data=ToxAngIgBF) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-BF,sum),y=BF,color = BF < 0, fill = BF < 0) +
  scale_color_manual(values=c("#38598CFF","grey40")) + scale_fill_manual(values=c("#38598CFF","grey40"))+ 
  labs(y="BF Hypoth vs Average", title="(Anguimorph,Iguania)")  
AI
SA <- ggplot(data=ToxSnAngBF) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-BF,sum),y=BF,color = BF < 0, fill = BF < 0) +
  scale_color_manual(values=c("#2BB07FFF","grey40")) + scale_fill_manual(values=c("#2BB07FFF","grey40"))+ 
  labs(y="BF Hypoth vs Average", title="(Snakes,Anguimorph)")  
SA
SI <- ggplot(data=ToxSnIgBF) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-BF,sum),y=BF,color = BF < 0, fill = BF < 0) +
  scale_color_manual(values=c("#C2DF23FF","grey40")) + scale_fill_manual(values=c("#C2DF23FF","grey40"))+ 
  labs(y="BF Hypoth vs Average", title="(Snakes,Iguania)")  
SI



all <- ggarrange(SC, AI, SA, SI, ncol=1, nrow=4, common.legend = TRUE, legend=FALSE)
all
ggsave("Streicher_BF_Classic_All.jpg", plot=all, width = 20, height = 40, units = "cm")



AI <- ggplot(data=ToxAngIgBFtx) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-BF,sum),y=BF,color = BF < 0, fill = BF < 0) +
  scale_color_manual(values=c("#38598CFF","grey40")) + scale_fill_manual(values=c("#38598CFF","grey40"))+ 
  labs(y="BF Hypoth vs Average", title="(Anguimorph,Iguania)")  
AI
SA <- ggplot(data=ToxSnAngBFtx) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-BF,sum),y=BF,color = BF < 0, fill = BF < 0) +
  scale_color_manual(values=c("#2BB07FFF","grey40")) + scale_fill_manual(values=c("#2BB07FFF","grey40"))+ 
  labs(y="BF Hypoth vs Average", title="(Snakes,Anguimorph)")  
SA
SI <- ggplot(data=ToxSnIgBFtx) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-BF,sum),y=BF,color = BF < 0, fill = BF < 0) +
  scale_color_manual(values=c("#C2DF23FF","grey40")) + scale_fill_manual(values=c("#C2DF23FF","grey40"))+ 
  labs(y="BF Hypoth vs Average", title="(Snakes,Iguania)")  
SI

all <- ggarrange(AI, SA, SI, ncol=1, nrow=3, common.legend = TRUE, legend=FALSE)
all
ggsave("Streicher_BF_Classic_Tox.jpg", plot=all, width = 20, height = 30, units = "cm")


# Not saved, histograms of locus types. should be split into different hypotheses for full information effect 

all <- ggplot(allBF, aes(x=BF, fill=LocusType, color=LocusType)) + 
  geom_histogram(binwidth = 5, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) + ylim(0,200)

all

all.trunc <- ggplot(allBF, aes(x=BF, fill=LocusType, color=LocusType)) + 
  geom_histogram(binwidth = 5, alpha=0.5, position="dodge") +
  geom_vline(xintercept=c(0,10),color=c("grey55", "black"), linetype="dashed", size=0.5) + 
  xlim(-200,200)

all
all.trunc 