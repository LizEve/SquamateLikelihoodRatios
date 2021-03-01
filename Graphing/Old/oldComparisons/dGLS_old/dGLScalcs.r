rm(list=ls())
library(viridis)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(gridExtra)
library(ggpubr)


# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
#scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

#scale_color_viridis(discrete=TRUE) 

dGLSofAvgToxNum <- function(numeratorConstraint,maxLikes){
  # Rename column name for numberator in calculation. 
  colnames(maxLikes)[colnames(maxLikes)==numeratorConstraint] <- "numeratorConstraint"
  
  # Start a list of dGLS values for all loci in dataset
  numLoci <- length(maxLikes$numeratorConstraint)
  x <- "ToxAll"
  # Make lists and add to data frame with locus names, constraint name, and 0 values
  numeratorList <- rep(c(x), times=numLoci)
  outdGLSList <- vector(mode="numeric",length=numLoci)
  
  # Make dataframe 
  outdGLSdf <- subset(maxLikes, select = c('Locus'))
  outdGLSdf$Hypothesis <- numeratorList
  outdGLSdf$dGLSavg <- outdGLSList
  outdGLSdf$dGLSmax <- outdGLSList
  
  # create a dataset of alternate constraints only to average over
  drops <- c("Locus","numeratorConstraint")
  df <- maxLikes[,!(names(maxLikes) %in% drops)]
  
  # for locus in dataset 
  for (i in 1:nrow(outdGLSdf)){
    
    # Get all alternate maxinal likehoods 
    othermaxLikes <- as.numeric(df[i,])
    
    # Maximum of alternate maxLH
    a <- max(othermaxLikes)
    
    # Calculate average across alternate constraints
    avgOthermaxLikes <- mean(othermaxLikes)
    # Calculate average across alternate constraints
    avgOtherMargLikes <- a+log(sum(exp(othermaxLikes-a)))-log(length(names(df)))
    
    # Calculate dGLS based on average of other LH values
    outdGLSdf[i,3] <- avgOthermaxLikes-maxLikes$numeratorConstraint[i]
    # Calculate dGLS based on maximum of other LH values
    outdGLSdf[i,4] <- a-maxLikes$numeratorConstraint[i]}
  
  return(outdGLSdf)
}

dGLSofAvg <- function(numeratorConstraint,maxLikes){
  
  # Rename column name for numberator in calculation. 
  colnames(maxLikes)[colnames(maxLikes)==numeratorConstraint] <- "numeratorConstraint"
  
  # Start a list of dGLS values for all loci in dataset
  numLoci <- length(maxLikes$numeratorConstraint)
  
  # Make lists and add to data frame with locus names, constraint name, and 0 values
  numeratorList <- rep(c(numeratorConstraint), times=numLoci)
  outdGLSList <- vector(mode="numeric",length=numLoci)
  
  # Make dataframe 
  outdGLSdf <- subset(maxLikes, select = c('Locus'))
  outdGLSdf$Hypothesis <- numeratorList
  outdGLSdf$dGLSavg <- outdGLSList
  outdGLSdf$dGLSmax <- outdGLSList
    
  # create a dataset of alternate constraints only to average over
  drops <- c("Locus","numeratorConstraint")
  df <- maxLikes[,!(names(maxLikes) %in% drops)]
  
  # for locus in dataset 
  for (i in 1:nrow(outdGLSdf)){
    
    # Get all alternate maxinal likehoods 
    othermaxLikes <- as.numeric(df[i,])
    
    # Maximum of alternate maxLH
    a <- max(othermaxLikes)
    
    # Calculate average across alternate constraints
    avgOthermaxLikes <- mean(othermaxLikes)
    # Calculate average across alternate constraints
    avgOtherMargLikes <- a+log(sum(exp(othermaxLikes-a)))-log(length(names(df)))
    
    # Calculate dGLS based on average of other LH values
    outdGLSdf[i,3] <- maxLikes$numeratorConstraint[i]-avgOthermaxLikes
   # Calculate dGLS based on maximum of other LH values
    outdGLSdf[i,4] <- maxLikes$numeratorConstraint[i]-a}
  
  return(outdGLSdf)
}

writeCSV <- function (outData,outFile){
  write.csv(outData, file=outFile, row.names=FALSE)
}


setwd("/Users/ChatNoir/Projects/Squam/Graphs/dGLS")
mLdGLS <- read.table("SquamNTGdGLSSummary.txt",header=TRUE)

# Read in table of marginal Likelihoods
mL <- read.table("SquamNTGdGLSSummary.txt",header=TRUE)
mLTox <- mL[,!(names(mL) %in% c("Sclero"))]
names(mL)
mLaiS <- mL[,!(names(mL) %in% c("ToxAngIg"))]
mLsaS <- mL[,!(names(mL) %in% c("ToxSnIg"))]
mLsiS <- mL[,!(names(mL) %in% c("ToxSnAng"))]


# Get list of dGLS for each relationship, using an average of all other constraints as a composite hypothesis. 

SclerodGLS <- dGLSofAvg("Sclero",mL)
SaidGLS <- dGLSofAvg("Sclero",mLaiS)
SsadGLS <- dGLSofAvg("Sclero",mLsaS)
SsidGLS <- dGLSofAvg("Sclero",mLsiS)
# differences only within tox constraints 
ToxAngIgdGLStx  <- dGLSofAvg("ToxAngIg",mLTox)
ToxSnAngdGLStx  <- dGLSofAvg("ToxSnAng",mLTox) 
ToxSnIgdGLStx <- dGLSofAvg("ToxSnIg",mLTox)
# all
ToxAngIgdGLS  <- dGLSofAvg("ToxAngIg",mL)
ToxSnAngdGLS  <- dGLSofAvg("ToxSnAng",mL) 
ToxSnIgdGLS <- dGLSofAvg("ToxSnIg",mL)

# Combine and factor 
alldGLS <- bind_rows(SclerodGLS,ToxAngIgdGLS,ToxSnAngdGLS,ToxSnIgdGLS)
toxdGLS <- bind_rows(ToxAngIgdGLStx,ToxSnAngdGLStx,ToxSnIgdGLStx)
alldGLS$Hypothesis <- factor(alldGLS$Hypothesis, levels=c("Sclero","ToxAngIg","ToxSnAng","ToxSnIg"))
toxdGLS$Hypothesis <- factor(toxdGLS$Hypothesis, levels=c("ToxAngIg","ToxSnAng","ToxSnIg"))
toxdGLS <- toxdGLS %>% separate(Locus, c("LocusType","LocusNumber"), "-", remove = FALSE)

# Absolute values, averages of differences 

mLablsAve <- mL %>% mutate(absAvg=(abs(Sclero-ToxAngIg)+
                        abs(Sclero-ToxSnAng)+
                        abs(Sclero-ToxSnIg)+
                        abs(ToxAngIg-ToxSnAng)+
                        abs(ToxAngIg-ToxSnIg)+
                        abs(ToxSnAng-ToxSnIg))/6) %>% select(Locus,absAvg)



# Support counts 


# SAve HERE 

#writeCSV(tvsSupport,"tvsdGLSavg_support_avg.csv")
#writeCSV(toxSupport,"toxdGLSavg_support_avg.csv")
#writeCSV(allSupport,"alldGLSavg_support.csv")

# Loci Lists 

# Loci sig support for Sclero over Tox 
# Loci sig support for Tox over Sclero 
# Loci sig support for AI over other Tox 
# Loci sig support for SA over other Tox
# Loci sig support for SI over other Tox

alldGLS.SC.loci <- alldGLS %>% select(Locus, dGLSavg, Hypothesis, Percentile, locusRank) %>% 
  filter(Hypothesis=="Sclero") %>% 
  filter(dGLSavg >= 0.5) %>% select(-Hypothesis)
alldGLS.AI.loci <- alldGLS %>% select(Locus, dGLSavg,Hypothesis,Percentile, locusRank) %>% 
  filter(Hypothesis=="ToxAngIg") %>% 
  filter(dGLSavg >= 0.5) %>% select(-Hypothesis) 
alldGLS.SA.loci <- alldGLS %>% select(Locus, dGLSavg,Hypothesis,Percentile, locusRank) %>% 
  filter(Hypothesis=="ToxSnAng") %>% 
  filter(dGLSavg >= 0.5) %>% select(-Hypothesis) 
alldGLS.SI.loci <- alldGLS %>% select(Locus, dGLSavg,Hypothesis,Percentile, locusRank) %>% 
  filter(Hypothesis=="ToxSnIg") %>% 
  filter(dGLSavg >= 0.5) %>% select(-Hypothesis)

toxdGLS.AI.loci <- toxdGLS %>% select(Locus, dGLSavg,Hypothesis,Percentile) %>% 
  filter(Hypothesis=="ToxAngIg") %>% 
  filter(dGLSavg >= 0.5) %>% select(-Hypothesis) 
toxdGLS.SA.loci <- toxdGLS %>% select(Locus, dGLSavg,Hypothesis,Percentile) %>% 
  filter(Hypothesis=="ToxSnAng") %>% 
  filter(dGLSavg >= 0.5) %>% select(-Hypothesis) 
toxdGLS.SI.loci <- toxdGLS %>% select(Locus, dGLSavg,Hypothesis,Percentile) %>% 
  filter(Hypothesis=="ToxSnIg") %>% 
  filter(dGLSavg >= 0.5) %>% select(-Hypothesis)

#writeCSV(tvsdGLS.SC.loci,"tvsdGLS.SC.loci.avg.csv")
#writeCSV(tvsdGLS.TX.loci,"tvsdGLS.TX.loci.avg.csv")
#writeCSV(toxdGLS.AI.loci,"tvsdGLS.AI.loci.avg.csv")
#writeCSV(toxdGLS.SA.loci,"tvsdGLS.SA.loci.avg.csv")
#writeCSV(toxdGLS.SI.loci,"tvsdGLS.SI.loci.avg.csv")

############# Bar graphs Pair dGLS
dGLSpairSupport <- pairdGLS %>% group_by(Hypothesis,Support) %>% dplyr::summarise(dGLSpairCounts=n())
dGLSpairSupport$Support <- factor(dGLSpairSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))


############## Plots ############## ############## ############## ############## ############## 


quartz()

max <- ggplot(data=ToxSnIgdGLStx) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 1) + 
  aes(x=reorder(Locus,-dGLSmax,sum),y=dGLSmax,color = dGLSmax < 0, fill = dGLSmax < 0) +
  scale_color_manual(values=c("#990000","#482173FF")) + 
  scale_fill_manual(values=c("#990000","#482173FF"))+ 
  labs(x="",y="Support Value",size=24) +
  theme_classic() + 
  theme(
    axis.line=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x= element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size=24, color="black"),
    text = element_text(size=30),
    legend.position = "none", 
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) 
max + scale_y_continuous(breaks=seq(20,20, 5))
ggsave("ToyExample_dGLSPlot2.pdf", plot=max, 
       width = 11, height = 4, units = "in", device = 'pdf',
       bg = "transparent")


############## Saved Plots ############## ############## ############## ############## ############## 
tvsList <- c("orange","#482173FF")
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
#show_col(myList)

# Singhal_dGLS_Histogram_TvS

myList <- c("orange","#482173FF")
avg <- ggplot(tvsdGLS, aes(x=dGLSavg, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  labs(x="dGLS Hypoth vs average")
max <- ggplot(tvsdGLS, aes(x=dGLSmax, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  labs(x="dGLS Hypoth vs maximum")
all <- ggarrange(avg,max, ncol=2, nrow=1, common.legend = TRUE, legend="right")
all
ggsave("Singhal_dGLS_Histogram_TvS_avg.jpg", plot=avg, width = 20, height = 20, units = "cm")

# Singhal_dGLS_Histogram_All 

myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
avg <- ggplot(alldGLS, aes(x=dGLSavg, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  labs(x="dGLS Hypoth vs average")
max <- ggplot(alldGLS, aes(x=dGLSmax, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  labs(x="dGLS Hypoth vs maximum")
all <- ggarrange(avg,max, ncol=2, nrow=1, common.legend = TRUE, legend="right")
all
ggsave("Singhal_dGLS_Histogram_All_avgvmax.jpg", plot=all, width = 30, height = 20, units = "cm")

all <- ggplot(alldGLS, aes(x=dGLSavg, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=myList)+scale_fill_manual(values=myList)
all
ggsave("Singhal_dGLS_Histo_all.jpg", plot=all, width = 40, height = 20, units = "cm")


# Singhal_dGLS_Violin_avg 

tvs <- ggplot(tvsdGLS, aes(x=Hypothesis, y=dGLSavg, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + scale_color_manual(values=tvsList) + stat_summary(fun.data=data_summary)
tox <- ggplot(toxdGLS, aes(x=Hypothesis, y=dGLSavg, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + scale_color_manual(values=myList) + stat_summary(fun.data=data_summary)
tox.trunc <- ggplot(toxdGLS, aes(x=Hypothesis, y=dGLSavg, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + scale_color_manual(values=myList) + stat_summary(fun.data=data_summary) +
  ylim(-2,2)
tox.trunc

ggsave("Singhal_dGLS_Violin_TvS_avg.jpg", plot=tvs, width = 20, height = 20, units = "cm")
ggsave("Singhal_dGLS_Violin_Tox_avg.jpg", plot=tox, width = 20, height = 20, units = "cm")
ggsave("Singhal_dGLS_Violin_Tox_avgtrunc.jpg", plot=tox.trunc, width = 20, height = 20, units = "cm")


# Singhal_dGLS_Classic_avg

tvsList <- c("orange","#482173FF")

max <- ggplot(data=SclerodGLS) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-dGLSmax,sum),y=dGLSmax,color = dGLSmax < 0, fill = dGLSmax < 0) +
  scale_color_manual(values=c("orange","#482173FF")) + scale_fill_manual(values=c("orange","#482173FF"))+ 
  labs(y="dGLS Hypoth vs maximum")
avg <- ggplot(data=SclerodGLS) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-dGLSavg,sum),y=dGLSavg,color = dGLSavg < 0, fill = dGLSavg < 0) +
  scale_color_manual(values=c("orange","#482173FF")) + scale_fill_manual(values=c("orange","#482173FF"))+
  labs(y="dGLS Hypoth vs average")
avg

all <- ggarrange(avg,max, ncol=1, nrow=2, common.legend = TRUE, legend=FALSE)
all
ggsave("Singhal_dGLS_Classic_TvS_avg.jpg", plot=avg, width = 30, height = 20, units = "cm")


toxList <- c("#38598CFF","#2BB07FFF","#C2DF23FF")

avgAI <- ggplot(data=ToxAngIgdGLStx) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-dGLSavg,sum),y=dGLSavg,color = dGLSavg < 0, fill = dGLSavg < 0) +
  scale_color_manual(values=c("#38598CFF","grey40")) + scale_fill_manual(values=c("#38598CFF","grey40"))+ 
  labs(y="dGLS Hypoth vs Average", title="(Anguimorph,Iguania)")  + scale_y_continuous(breaks=seq(-17,17, 2), limits=c(-17,17)) 
#avgAI
avgSA <- ggplot(data=ToxSnAngdGLStx) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-dGLSavg,sum),y=dGLSavg,color = dGLSavg < 0, fill = dGLSavg < 0) +
  scale_color_manual(values=c("#2BB07FFF","grey40")) + scale_fill_manual(values=c("#2BB07FFF","grey40"))+ 
  labs(y="dGLS Hypoth vs Average", title="(Snakes,Anguimorph)")  + scale_y_continuous(breaks=seq(-17,17, 2), limits=c(-17,17))
#maxSA
avgSI <- ggplot(data=ToxSnIgdGLStx) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  aes(x=reorder(Locus,-dGLSavg,sum),y=dGLSavg,color = dGLSavg < 0, fill = dGLSavg < 0) +
  scale_color_manual(values=c("#C2DF23FF","grey40")) + scale_fill_manual(values=c("#C2DF23FF","grey40"))+ 
  labs(y="dGLS Hypoth vs Average", title="(Snakes,Iguania)")  + scale_y_continuous(breaks=seq(-17,17, 2), limits=c(-17,17))
#avgSI


all <- ggarrange(avgAI, avgSA, avgSI, ncol=1, nrow=3, common.legend = TRUE, legend=FALSE)
all
ggsave("Singhal_dGLS_Classic_Tox_avg.jpg", plot=all, width = 20, height = 40, units = "cm")



