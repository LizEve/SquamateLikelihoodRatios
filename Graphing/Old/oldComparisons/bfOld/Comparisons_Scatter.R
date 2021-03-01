rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("/Users/ChatNoir/Projects/Squam/Graphs")


load("/Users/ChatNoir/Projects/Squam/Graphs/BF/Singhal_BFdata.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/Singhal_dGLSdata.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/GGI/Singhal_GGIdata.RData")

#load("/Users/ChatNoir/Projects/Squam/Graphs/BF/Streicher_BFdata.RData")
#load("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/Streicher_dGLSdata.RData")
#load("/Users/ChatNoir/Projects/Squam/Graphs/GGI/Streicher_GGIdata.RData")

allBF <- allBF %>% select(-LocusNumber,-locusRank,-Percentile)%>%
  group_by(Hypothesis) %>% 
  mutate(Percentile = ntile(desc(BF),100)) %>% 
  mutate(Decile = ntile(desc(BF),10)) %>% 
  mutate(locusRank=dense_rank(desc(BF))) %>% 
  rename(BF.Percentile=Percentile,BF.Decile=Decile, 
         BF.locusRank = locusRank)
# This will give the same rank number if the value is identicle in two cells 

alldGLS <- alldGLS %>%  select(-dGLSmax,-locusRank,-Percentile) %>%
  group_by(Hypothesis) %>% 
  mutate(Percentile = ntile(desc(dGLSavg),100)) %>% 
  mutate(Decile = ntile(desc(dGLSavg),10)) %>% 
  mutate(locusRank=dense_rank(desc(dGLSavg))) %>% 
  rename(dGLS.Percentile=Percentile, dGLS.Decile=Decile, 
         dGLS.locusRank = locusRank) %>%
  mutate(Support = str_replace(Support, " < -0.5", "")) %>% 
  mutate(Support = str_replace(Support, " > 0.5", "")) 

allGGI <- allAU %>% rename(Locus=LOCUS) %>% select(-item)%>%
  group_by(Hypothesis) %>% 
  mutate(Percentile = ntile(desc(au),100)) %>% 
  mutate(Decile = ntile(desc(au),10)) %>% 
  mutate(locusRank=dense_rank(desc(au))) %>% 
  rename(GGI.Percentile=Percentile, GGI.Decile=Decile, 
         GGI.locusRank = locusRank, GGI.Rank=rank, Support=Support, SupportAU=SupportAU)


######## Comparison #############
#dGLS vs BF
#BF vs GGI
#dGLS vs GGI 


all <- merge(merge(allBF, alldGLS, by = c("Locus","Hypothesis")),allGGI, by = c("Locus","Hypothesis"))

SC <- filter(all,Hypothesis == "Sclero")
AI <- filter(all,Hypothesis == "ToxAngIg")
SA <- filter(all,Hypothesis == "ToxSnAng")
SI <- filter(all,Hypothesis == "ToxSnIg")

quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")

cc <- c("#38598CFF")
# color by locus rank 
z <- ggplot(AI, aes(x=BF,y=dGLSavg)) + 
  geom_point(alpha=0.5, color=cc, size=3) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(size=20, color="black"),
    axis.text.y = element_text(size=20, color="black"),
    text = element_text(size=30),
    legend.position = "none", 
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.key=element_blank()# get rid of legend panel bg
  ) 
z
ggsave("Singhal_Scatter_AI.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")

# color by GGI rank 
z <- ggplot(AI, aes(x=BF,y=dGLSavg)) + 
  geom_point(size=3,alpha=0.5, aes(color=factor(GGI.Rank))) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=24, color="black"),
    axis.text.y = element_text(size=24, color="black"),
    text = element_text(size=30),
    legend.position = "top", 
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.key=element_blank()
    ) + 
  labs(color = "GGI rank") 
z
ggsave("Singhal_Scatter_AI_byGGIrank.pdf", plot=z, 
       width = 7, height = 6,
       units = "in", device = 'pdf',
       bg = "transparent")


# color by pval 
z <- ggplot(AI, aes(x=BF,y=dGLSavg)) + 
  geom_point(alpha=0.5, aes(color=factor(au))) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(size=24, color="black"),
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
z
ggsave("Singhal_Scatter_AI_pval.pdf", plot=z, 
       width = 6, height = 6,
       units = "in", device = 'pdf',
       bg = "transparent")


z <- ggplot(all, aes(x=BF,y=au)) + 
  geom_point(alpha=0.5, aes(color=factor(GGI.Rank))) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5)
z
ggsave("Singhal_Scatter_all_byGGIrank.jpg", plot=z, width = 10, height = 10, units = "cm")

cc="grey29"
a <- ggplot(all, aes(x=BF,y=dGLSavg)) + 
  geom_point(alpha=0.5,color=cc) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5)
b <- ggplot(all, aes(x=BF,y=au)) + 
  geom_point(alpha=0.5, aes(color=factor(GGI.Rank))) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.95),color=c("black"), linetype="dashed", size=0.5)
c <- ggplot(all, aes(x=au,y=dGLSavg)) + 
  geom_point(alpha=0.5, aes(color=factor(GGI.Rank))) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(0.95),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5)
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="bottom", align="v",common.legend = TRUE)
x
#ggsave("Streicher_Scatter_all.jpg", plot=x, width = 30, height = 10, units = "cm")


a <- ggplot(all, aes(x=BF,y=dGLSavg)) + 
  geom_point(alpha=0.5,aes(color=all$Hypothesis)) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  scale_color_manual(values=myList)
b <- ggplot(all, aes(x=BF,y=au)) + 
  geom_point(alpha=0.5, aes(color=all$Hypothesis)) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.95),color=c("black"), linetype="dashed", size=0.5) +
  scale_color_manual(values=myList)
c <- ggplot(all, aes(x=au,y=dGLSavg)) + 
  geom_point(alpha=0.5, aes(color=all$Hypothesis)) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(0.95),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  scale_color_manual(values=myList)
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="bottom", align="v",common.legend = TRUE)
x
ggsave("Streicher_Scatter_allbyHypoth.jpg", plot=a, width = 30, height = 10, units = "cm")

a <- ggplot(SC, aes(x=BF,y=dGLSavg)) + 
  geom_point(alpha=0.5,color=myList[1]) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5)
b <- ggplot(SC, aes(x=BF,y=au)) + 
  geom_point(alpha=0.5, aes(color=factor(GGI.Rank))) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.95),color=c("black"), linetype="dashed", size=0.5)
c <- ggplot(SC, aes(x=au,y=dGLSavg)) + 
  geom_point(alpha=0.5, aes(color=factor(GGI.Rank))) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(0.95),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5)
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="bottom", align="v",common.legend = TRUE)
x
ggsave("Streicher_Scatter_SC.jpg", plot=x, width = 30, height = 10, units = "cm")



# Color by locus type
a <- ggplot(all, aes(x=BF.locusRank,y=dGLS.locusRank,color=LocusType)) + 
  geom_point(alpha=0.5) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.3) 
b <- ggplot(all, aes(x=BF.locusRank,y=GGI.locusRank,color=LocusType)) + 
  geom_point(alpha=0.5) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.3) 
c <- ggplot(all, aes(x=GGI.locusRank,y=dGLS.locusRank,color=LocusType)) + 
  geom_point(alpha=0.5) +
  scale_y_reverse()+scale_x_reverse()+ theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.3) 
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="bottom", align="v")
x
ggsave("Singhal_LocusRank_all_LocusType.jpg", plot=x, width = 30, height = 10, units = "cm")

# Color by hypothesis 
cc <- myList[1]
a <- ggplot(SC, aes(x=BF.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
b <- ggplot(SC, aes(x=BF.locusRank,y=GGI.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
c <- ggplot(SC, aes(x=GGI.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse()+ theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="none", align="v")
x
#ggsave("Streicher_LocusRank_SC.jpg", plot=x, width = 30, height = 10, units = "cm")
ggsave("Singhal_LocusRank_SC.jpg", plot=x, width = 30, height = 10, units = "cm")


cc <- myList[2]
a <- ggplot(AI, aes(x=BF.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
b <- ggplot(AI, aes(x=BF.locusRank,y=GGI.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
c <- ggplot(AI, aes(x=GGI.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse()+ theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="none", align="v")
x
#ggsave("Streicher_LocusRank_AI.jpg", plot=x, width = 30, height = 10, units = "cm")
ggsave("Singhal_LocusRank_AI.jpg", plot=x, width = 30, height = 10, units = "cm")

cc <- myList[3]
a <- ggplot(SA, aes(x=BF.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
b <- ggplot(SA, aes(x=BF.locusRank,y=GGI.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
c <- ggplot(SA, aes(x=GGI.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse()+ theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="none", align="v")
x

#ggsave("Streicher_LocusRank_SA.jpg", plot=x, width = 30, height = 10, units = "cm")
ggsave("Singhal_LocusRank_SA.jpg", plot=x, width = 30, height = 10, units = "cm")


cc <- myList[4]
a <- ggplot(SI, aes(x=BF.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
b <- ggplot(SI, aes(x=BF.locusRank,y=GGI.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
c <- ggplot(SI, aes(x=GGI.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse()+ theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="none", align="v")
x

#ggsave("Streicher_LocusRank_SI.jpg", plot=x, width = 30, height = 10, units = "cm")
ggsave("Singhal_LocusRank_SI.jpg", plot=x, width = 30, height = 10, units = "cm")





##### Counts #########



# Count strong support loci 
BFallSupport <- allBF %>% group_by(Hypothesis,BF.Support) %>% dplyr::summarise(BFallCounts=n())
BFallSupport$BFallSupport <- factor(BFallSupport$BF.Support, levels=c("Strong_Against","Ambiguous","Strong"))

GGIallSupport <- allGGI %>% group_by(Hypothesis,GGI.Support) %>% dplyr::summarise(GGIallCounts=n())
GGIallSupport$GGIallSupport <- factor(GGIallSupport$GGI.Support, levels=c("Strong_Against","Ambiguous","Strong"))

dGLSallSupport <- alldGLS %>% group_by(Hypothesis,dGLS.Support) %>% dplyr::summarise(dGLSallCounts=n())
dGLSallSupport$dGLSallSupportavg <- factor(dGLSallSupport$dGLS.Support, levels=c("Strong_Against","Ambiguous","Strong"))

# Count shared loci with support by hypothesis 

x <- subset(allAU,select=c("LOCUS")) %>% rename(Locus=LOCUS) %>% distinct()
b <- allBF %>% filter(Hypothesis == "Sclero") 
d <- alldGLS %>% filter(Hypothesis == "Sclero") 
g <- allGGI %>% filter(Hypothesis == "Sclero")
SC <- left_join(left_join(left_join(x,b, by = c("Locus")),d, by = c("Locus")),g,by=c("Locus"))

x <- subset(allAU,select=c("LOCUS")) %>% rename(Locus=LOCUS) %>% distinct()
b <- allBF %>% filter(Hypothesis == "Sclero") %>% filter(BF.Support=="Strong")
d <- alldGLS %>% filter(Hypothesis == "Sclero") %>% filter(dGLS.Support=="Strong")
g <- allGGI %>% filter(Hypothesis == "Sclero")%>% filter(GGI.Support=="Strong")
SC.strong <- left_join(left_join(left_join(x,b, by = c("Locus")),d, by = c("Locus")),g,by=c("Locus"))



SC.all.loci <- SC %>% 
  mutate(methodType=case_when(BF >= 10 & dGLSavg >= 0.5 & GGI.Rank == 1 ~ "all",
                              BF >= 10 & dGLSavg >= 0.5 & is.na(au) ~ "BF,dGLS",
                              BF >= 10 & is.na(dGLSavg) & GGI.Rank == 1 ~ "BF,GGI",
                              is.na(BF) & dGLSavg >= 0.5 & GGI.Rank == 1 ~ "dGLS,GGI",
                              BF >= 10  ~ "onlyBF", 
                              dGLSavg >= 0.5 ~ "onlydGLS",
                              au >= 0 ~ "onlyGGI",
                              TRUE ~ "no.strong.support"))

SC.all.loci.sum <- SC.all.loci %>%   group_by(methodType) %>% summarize(loci=n()) %>% filter(!(methodType=='no.strong.support'))







#### Graph 
quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
sc.all <- ggplot(SC.all.loci.sum, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[1]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") + theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250) 

ai.all <- ggplot(AI.all.loci.sum, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[2]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250)

sa.all <- ggplot(SA.all.loci.sum, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[3]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250)

sa.each <- ggplot(SA.each.loci, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[3]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250)

sa <- ggarrange(sa.each,sa.all, ncol=2, nrow=1, common.legend = TRUE, legend="right")
sa

si.all <- ggplot(SI.all.loci, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[4]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250)

si.each <- ggplot(SI.each.loci, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[4]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250)

si <- ggarrange(si.each,si.all, ncol=2, nrow=1, common.legend = TRUE, legend="right")
si

s <- ggarrange(sc,ai,sa,si, ncol=1, nrow=4)
s
ggsave("Singhal_LociMethods.jpg", plot=s, width = 40, height = 30, units = "cm")





#### VENNN

source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
library(tidyverse)
library(ggforce)

# Fix later 

df.vdc <- as.data.frame(SC.all.loci.sum) %>%
  mutate(x= c(0,0.8,),
         y=(0,0.5,))
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))

df.venn <- data.frame(x = c(0, 0.866, -0.866),
                      y = c(1, -0.5, -0.5),
                      labels = c('BF', 'dGLS', 'GGI'))
ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom') +
  scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
  scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$loci, size = 5)

# More counts 
