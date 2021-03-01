rm(list=ls())
library(viridis)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(gridExtra)
library(ggpubr)

# To use for fills, add
#scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

#scale_color_viridis(discrete=TRUE) 

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
mLdGLS <- read.table("StreicherdGLSSummary.txt",header=TRUE)

# Read in table of marginal Likelihoods
mL <- read.table("StreicherdGLSSummary.txt",header=TRUE)
mLTox <- mL[,!(names(mL) %in% c("Sclero"))]
names(mL)
mLaiS <- mL[,!(names(mL) %in% c("ToxAngIg"))]
mLsaS <- mL[,!(names(mL) %in% c("ToxSnIg"))]
mLsiS <- mL[,!(names(mL) %in% c("ToxSnAng"))]



# Get list of dGLS for each relationship, using an average of all other constraints as a composite hypothesis. 

SclerodGLS <- dGLSofAvg("Sclero",mL)
SaiBF <- dGLSofAvg("Sclero",mLaiS)
SsaBF <- dGLSofAvg("Sclero",mLsaS)
SsiBF <- dGLSofAvg("Sclero",mLsiS)
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

toxdGLS <- mutate(toxdGLS, Supportavg = case_when(dGLSavg >= 0.5 ~ "Strong > 0.5",
                                                  dGLSavg > -0.5 & dGLSavg < 0.5 ~ "Ambiguous",
                                                  dGLSavg <= -0.5 ~ "Strong_Against < -0.5")) %>%
  mutate(Percentile = ntile(desc(dGLSavg),100))%>% mutate(locusRank=dense_rank(desc(dGLSavg)))
alldGLS <- mutate(alldGLS, Supportavg = case_when(dGLSavg >= 0.5 ~ "Strong > 0.5",
                                                  dGLSavg > -0.5 & dGLSavg < 0.5 ~ "Ambiguous",
                                                  dGLSavg <= -0.5 ~ "Strong_Against < -0.5")) %>%
  mutate(Percentile = ntile(desc(dGLSavg),100))%>% mutate(locusRank=dense_rank(desc(dGLSavg)))



toxSupport <- toxdGLS %>% group_by(Hypothesis,Supportavg) %>% dplyr::summarise(counts=n())
toxSupport$Supportavg <- factor(toxSupport$Supportavg, levels=c("Strong_Against < -0.5","Ambiguous","Strong > 0.5"))
allSupport <- alldGLS %>% group_by(Hypothesis,Supportavg) %>% dplyr::summarise(counts=n())
allSupport$Supportavg <- factor(allSupport$Supportavg, levels=c("Strong_Against < -0.5","Ambiguous","Strong > 0.5"))

# SAve R data here 



#writeCSV(tvsSupport,"Streicher_tvsdGLSavg_support_avg.csv")
#writeCSV(toxSupport,"Streicher_toxdGLSavg_support_avg.csv")
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

#writeCSV(tvsdGLS.SC.loci,"Streicher_tvsdGLS.SC.loci.avg.csv")
#writeCSV(tvsdGLS.TX.loci,"Streicher_tvsdGLS.TX.loci.avg.csv")
#writeCSV(toxdGLS.AI.loci,"Streicher_tvsdGLS.AI.loci.avg.csv")
#writeCSV(toxdGLS.SA.loci,"Streicher_tvsdGLS.SA.loci.avg.csv")
#writeCSV(toxdGLS.SI.loci,"Streicher_tvsdGLS.SI.loci.avg.csv")





############## Plots ############## ############## ############## ############## ############## 

quartz()







############## Saved Plots ############## ############## ############## ############## ############## 
tvsList <- c("orange","#482173FF")
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
#show_col(myList)

# Streicher_dGLS_Histogram_All

myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
avg <- ggplot(alldGLS, aes(x=dGLSavg, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  labs(x="dGLS Hypoth vs average")
avg.trunc <- ggplot(alldGLS, aes(x=dGLSavg, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  labs(x="dGLS Hypoth vs avg") + xlim(-15,15)
avg.trunc
all <- ggarrange(avg,avg.trunc, ncol=1, nrow=2, common.legend = TRUE, legend="right")
avg
quartz()
all
ggsave("Streicher_dGLS_Histogram_All_avg.jpg", plot=all, width = 20, height = 20, units = "cm")


# Streicher_dGLS_Histogram_Tox - not useful visual 

myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
txList <- c("#38598CFF","#2BB07FFF","#C2DF23FF")

avg <- ggplot(toxdGLS, aes(x=dGLSavg, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=txList) + scale_fill_manual(values=txList) + 
  labs(x="dGLS Hypoth vs average") 
max <- ggplot(toxdGLS, aes(x=dGLSmax, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=txList) + scale_fill_manual(values=txList) + 
  labs(x="dGLS Hypoth vs maximum") 
all <- ggarrange(avg,max, ncol=2, nrow=1, common.legend = TRUE, legend="right")
all
avg.trunc <- ggplot(toxdGLS, aes(x=dGLSavg, fill=Hypothesis, color=Hypothesis)) + 
  geom_histogram(binwidth = 0.5, alpha=1, position="dodge") +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_manual(values=txList) + scale_fill_manual(values=txList) + 
  labs(x="dGLS Hypoth vs average")  + xlim(-15,15)
avg.trunc
ggsave("Streicher_dGLS_Histogram_Tox_trunc.jpg", plot=avg.trunc, width = 30, height = 20, units = "cm")


# Streicher_dGLS_Violin_avg 
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
tvs <- ggplot(tvsdGLS, aes(x=Hypothesis, y=dGLSavg, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + scale_color_manual(values=tvsList) + stat_summary(fun.data=data_summary)
tox <- ggplot(toxdGLS, aes(x=Hypothesis, y=dGLSavg, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + scale_color_manual(values=myList) + stat_summary(fun.data=data_summary)
tox.trunc <- ggplot(toxdGLS, aes(x=Hypothesis, y=dGLSavg, color=Hypothesis)) + 
  geom_violin(trim=TRUE) + scale_color_manual(values=myList) + stat_summary(fun.data=data_summary) +
  ylim(-2,2)
tox.trunc

#ggsave("Streicher_dGLS_Violin_TvS_avg.jpg", plot=tvs, width = 20, height = 20, units = "cm")
#ggsave("Streicher_dGLS_Violin_Tox_avg.jpg", plot=tox, width = 20, height = 20, units = "cm")
#ggsave("Streicher_dGLS_Violin_Tox_avgtrunc.jpg", plot=tox.trunc, width = 20, height = 20, units = "cm")


# Streicher_dGLS_Classic_avg

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
ggsave("Streicher_dGLS_Classic_TvS_avg.jpg", plot=avg, width = 30, height = 20, units = "cm")


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
#ggsave("Streicher_dGLS_Classic_Tox_avg.jpg", plot=all, width = 20, height = 40, units = "cm")



