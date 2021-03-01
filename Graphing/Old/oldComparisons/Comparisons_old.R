rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("/Users/ChatNoir/Projects/Squam/Graphs")

load("/Users/ChatNoir/Projects/Squam/Graphs/BF/Singhal_BFdata.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/dGLSdata.avg.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/GGI/GGIdata.RData")

######## Loci in common #############
ls(pattern="loci",all.names=TRUE) 
ls(pattern="all",all.names=TRUE)

# number of loci with strong support for each hypothesis all, and within tox 
# are they the same loci?
# number of loci with absolute strong signal for dgls and bf 
#Loci with BF >10 supporting X hypoth 

allBF.SC.loci
allBF.AI.loci
allBF.SA.loci
allBF.SI.loci

alldGLS.SC.loci
alldGLS.AI.loci
alldGLS.SA.loci
alldGLS.SI.loci

allAU.r1.SC.loci <- allAU.r1.loci %>% filter(Hypothesis == "Sclero") %>% rename(Locus=LOCUS)
allAU.r1.AI.loci <- allAU.r1.loci %>% filter(Hypothesis == "ToxAngIg") %>% rename(Locus=LOCUS)
allAU.r1.SA.loci <- allAU.r1.loci %>% filter(Hypothesis == "ToxSnAng") %>% rename(Locus=LOCUS)
allAU.r1.SI.loci <- allAU.r1.loci %>% filter(Hypothesis == "ToxSnIg") %>% rename(Locus=LOCUS)
allAU.r1.95.loci

# Grab all loci. combine lists with strong support loci 
x <- subset(allAU,select=c("LOCUS")) %>% rename(Locus=LOCUS) %>% distinct()
SC <- left_join(x,allBF.SC.loci) %>% rename(BF.Percentile=Percentile) %>%
  left_join(alldGLS.SC.loci) %>% rename(dGLS.Percentile=Percentile) %>%
  left_join(allAU.r1.SC.loci) 

SC.each.loci <- SC %>% 
  mutate(methodType=case_when(BF >= 10  ~ "totalBF", 
                            dGLSavg >= 0.5 ~ "totaldGLS",
                            au >= 0 ~ "totalGGI",
                            is.na(BF) & is.na(dGLSavg) & is.na(au) ~ "no.strong.support" )) %>%
  group_by(methodType) %>% summarize(loci=n())

SC.all.loci <- SC %>% 
  mutate(methodType=case_when(BF >= 10 & dGLSavg >= 0.5 & au >= 0 ~ "all",
                                                     BF >= 10 & dGLSavg >= 0.5 & is.na(au) ~ "BF,dGLS",
                                                     BF >= 10 & is.na(dGLSavg) & au >= 0 ~ "BF,GGI",
                                                     is.na(BF) & dGLSavg >= 0.5 & au >= 0 ~ "dGLS,GGI",
                                                     BF >= 10  ~ "onlyBF", 
                                                     dGLSavg >= 0.5 ~ "onlydGLS",
                                                     au >= 0 ~ "onlyGGI",
                                                     TRUE ~ "no.strong.support")) %>%
  group_by(methodType) %>% summarize(loci=n()) %>% filter(methodType != "no.strong.support")

### AI
AI <- left_join(x,allBF.AI.loci) %>% rename(BF.Percentile=Percentile) %>%
  left_join(alldGLS.AI.loci) %>% rename(dGLS.Percentile=Percentile) %>%
  left_join(allAU.r1.AI.loci) 

AI.each.loci <- AI %>% 
  mutate(methodType=case_when(BF >= 10  ~ "totalBF", 
                              dGLSavg >= 0.5 ~ "totaldGLS",
                              au >= 0 ~ "totalGGI",
                              is.na(BF) & is.na(dGLSavg) & is.na(au) ~ "no.strong.support" )) %>%
  group_by(methodType) %>% summarize(loci=n())

AI.all.loci <- AI %>% 
  mutate(methodType=case_when(BF >= 10 & dGLSavg >= 0.5 & au >= 0 ~ "all",
                              BF >= 10 & dGLSavg >= 0.5 & is.na(au) ~ "BF,dGLS",
                              BF >= 10 & is.na(dGLSavg) & au >= 0 ~ "BF,GGI",
                              is.na(BF) & dGLSavg >= 0.5 & au >= 0 ~ "dGLS,GGI",
                              BF >= 10  ~ "onlyBF", 
                              dGLSavg >= 0.5 ~ "onlydGLS",
                              au >= 0 ~ "onlyGGI",
                              TRUE ~ "no.strong.support")) %>%
  group_by(methodType) %>% summarize(loci=n()) %>% filter(methodType != "no.strong.support")

### SA 
SA <- left_join(x,allBF.SA.loci) %>% rename(BF.Percentile=Percentile) %>%
  left_join(alldGLS.SA.loci) %>% rename(dGLS.Percentile=Percentile) %>%
  left_join(allAU.r1.SA.loci) 

SA.each.loci <- SA %>% 
  mutate(methodType=case_when(BF >= 10  ~ "totalBF", 
                              dGLSavg >= 0.5 ~ "totaldGLS",
                              au >= 0 ~ "totalGGI",
                              is.na(BF) & is.na(dGLSavg) & is.na(au) ~ "no.strong.support" )) %>%
  group_by(methodType) %>% summarize(loci=n())

SA.all.loci <- SA %>% 
  mutate(methodType=case_when(BF >= 10 & dGLSavg >= 0.5 & au >= 0 ~ "all",
                              BF >= 10 & dGLSavg >= 0.5 & is.na(au) ~ "BF,dGLS",
                              BF >= 10 & is.na(dGLSavg) & au >= 0 ~ "BF,GGI",
                              is.na(BF) & dGLSavg >= 0.5 & au >= 0 ~ "dGLS,GGI",
                              BF >= 10  ~ "onlyBF", 
                              dGLSavg >= 0.5 ~ "onlydGLS",
                              au >= 0 ~ "onlyGGI",
                              TRUE ~ "no.strong.support")) %>%
  group_by(methodType) %>% summarize(loci=n()) %>% filter(methodType != "no.strong.support")

### SI 

SI <- left_join(x,allBF.SI.loci) %>% rename(BF.Percentile=Percentile) %>%
  left_join(alldGLS.SI.loci) %>% rename(dGLS.Percentile=Percentile) %>%
  left_join(allAU.r1.SI.loci) 

SI.each.loci <- SI %>% 
  mutate(methodType=case_when(BF >= 10  ~ "totalBF", 
                              dGLSavg >= 0.5 ~ "totaldGLS",
                              au >= 0 ~ "totalGGI",
                              is.na(BF) & is.na(dGLSavg) & is.na(au) ~ "no.strong.support" )) %>%
  group_by(methodType) %>% summarize(loci=n())

SI.all.loci <- SI %>% 
  mutate(methodType=case_when(BF >= 10 & dGLSavg >= 0.5 & au >= 0 ~ "all",
                              BF >= 10 & dGLSavg >= 0.5 & is.na(au) ~ "BF,dGLS",
                              BF >= 10 & is.na(dGLSavg) & au >= 0 ~ "BF,GGI",
                              is.na(BF) & dGLSavg >= 0.5 & au >= 0 ~ "dGLS,GGI",
                              BF >= 10  ~ "onlyBF", 
                              dGLSavg >= 0.5 ~ "onlydGLS",
                              au >= 0 ~ "onlyGGI",
                              TRUE ~ "no.strong.support")) %>%
  group_by(methodType) %>% summarize(loci=n()) %>% filter(methodType != "no.strong.support")

#### Graph 
quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
sc.all <- ggplot(SC.all.loci, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[1]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") + theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250) 

sc.each <- ggplot(SC.each.loci, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[1]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") +  theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250)  

sc <- ggarrange(sc.each,sc.all, ncol=2, nrow=1, common.legend = TRUE, legend="right")
sc

ai.all <- ggplot(AI.all.loci, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[2]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250)

ai.each <- ggplot(AI.each.loci, aes(x=methodType,y=loci)) + 
  geom_bar(stat = "identity", fill=myList[2]) +
  geom_text(aes(label = loci), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,2250)

ai <- ggarrange(ai.each,ai.all, ncol=2, nrow=1, common.legend = TRUE, legend="right")
ai

sa.all <- ggplot(SA.all.loci, aes(x=methodType,y=loci)) + 
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

##### Counts #########

# Loci that support tox 283
GGIallSclero0.05 <- rall %>% group_by(allSclero0.05) %>% summarise(counts=n())
# 253
GGIaiSclero0.05 <- raiS %>% group_by(aiSSclero0.05) %>% summarise(counts=n()) 
# 268
GGIsaSclero0.05 <- rsaS %>% group_by(saSSclero0.05) %>% summarise(counts=n()) 
# 253
GGIsiSclero0.05 <- rsiS %>% group_by(siSSclero0.05) %>% summarise(counts=n()) 

# Count rank 1, and rank 1 with all others <= 0.05
GGIallRank1 <- allAU0.05 %>% group_by(Hypothesis,rank) %>% 
  filter(rank == 1) %>%  summarise(GGIallSupportR1=n()) %>% 
  mutate(rank = str_replace(rank, "1", "Strong")) %>% rename(Support=rank)
GGIallRank1.05 <- allAU0.05 %>% group_by(Hypothesis,rank) %>% 
  filter(rank == 1) %>% filter(lt0.05 == 3) %>% summarise(GGIallSupportR105=n()) %>% 
  mutate(rank = str_replace(rank, "1", "Strong"))

GGItoxRank1 <- toxAU0.05 %>% group_by(Hypothesis,rank) %>% 
  filter(rank == 1) %>%  summarise(GGItoxSupportR1=n()) %>% 
  mutate(rank = str_replace(rank, "1", "Strong")) %>% rename(Support=rank)
GGItoxRank1.05 <- toxAU0.05 %>% group_by(Hypothesis,rank) %>% 
  filter(rank == 1) %>% filter(lt0.05 == 3) %>% summarise(GGItoxSupportR105=n()) %>% 
  mutate(rank = str_replace(rank, "1", "Strong"))


# Count strong support loci 
BFallSupport <- allBF %>% group_by(Hypothesis,Support) %>% dplyr::summarise(BFallCounts=n())
BFallSupport$BFallSupport <- factor(BFallSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))
BFtoxSupport <- toxBF %>% group_by(Hypothesis,Support) %>% dplyr::summarise(BFtoxCounts=n())
BFtoxSupport$BFtoxSupport <- factor(BFtoxSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))

GGIallSupport <- allAU %>% group_by(Hypothesis,Support) %>% dplyr::summarise(GGIallCounts=n())
GGIallSupport$GGIallSupport <- factor(GGIallSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))
GGItoxSupport <- toxAU %>% group_by(Hypothesis,Support) %>% dplyr::summarise(GGItoxCounts=n())
GGItoxSupport$GGItoxSupport <- factor(GGItoxSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))


dGLSallSupport <- alldGLS %>% group_by(Hypothesis,Supportavg) %>% dplyr::summarise(dGLSallCounts=n()) %>% 
  mutate(Supportavg = str_replace(Supportavg, " < -0.5", "")) %>% 
  mutate(Supportavg = str_replace(Supportavg, " > 0.5", "")) %>% rename(Support=Supportavg)
dGLSallSupport$dGLSallSupportavg <- factor(dGLSallSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))
dGLStoxSupport <- toxdGLS %>% group_by(Hypothesis,Supportavg) %>% dplyr::summarise(dGLStoxCounts=n()) %>% 
  mutate(Supportavg = str_replace(Supportavg, " < -0.5", "")) %>% 
  mutate(Supportavg = str_replace(Supportavg, " > 0.5", "")) %>% rename(Support=Supportavg)
dGLStoxSupport$dGLStoxSupportavg <- factor(dGLStoxSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))

allSupport <- merge(GGIallSupport,merge(BFallSupport,dGLSallSupport,by=c("Hypothesis","Support")),by=c("Hypothesis","Support"))
allSupport <- left_join(allSupport, GGIallRank1,by=c("Hypothesis","Support"))
allSupport$Support <- factor(allSupport$Support, levels=c("Strong","Ambiguous","Strong_Against"))

toxSupport <- merge(GGItoxSupport,merge(BFtoxSupport,dGLStoxSupport,by=c("Hypothesis","Support")),by=c("Hypothesis","Support"))
toxSupport <- left_join(toxSupport, GGItoxRank1,by=c("Hypothesis","Support"))
toxSupport$Support <- factor(toxSupport$Support, levels=c("Strong","Ambiguous","Strong_Against"))



####### Bar_Singhall_Support Counts ###########

#quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
bf <- ggplot(allSupport, aes(x=Support,y=BFallCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = BFallCounts, group = Hypothesis), position = position_dodge(0.9),vjust = -0.3, size = 2) +
  labs(x="") + scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,3650)
dgls <- ggplot(allSupport, aes(x=Support,y=dGLSallCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = dGLSallCounts, group = Hypothesis), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") + scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,3650)
ggi <- ggplot(allSupport, aes(x=Support,y=GGIallCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGIallCounts, group = Hypothesis), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") + scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,3650)
ggir1 <- ggplot(allSupport, aes(x=Support,y=GGIallSupportR1)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGIallSupportR1, group = Hypothesis), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") + scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,3650)

a95 <- ggarrange(bf,dgls,ggi, ncol=1, nrow=3, common.legend = TRUE, legend="right") 
ar1 <- ggarrange(bf,dgls,ggir1, ncol=1, nrow=3, common.legend = TRUE, legend="right")

#ggsave("Singhal_SupportBarR1.jpg", plot=ar1, width = 20, height = 20, units = "cm")
#ggsave("Singhal_SupportBar95.jpg", plot=a95, width = 20, height = 30, units = "cm")

myList <- c("#38598CFF","#2BB07FFF","#C2DF23FF")
bf <- ggplot(toxSupport, aes(x=Support,y=BFtoxCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = BFtoxCounts, group = Hypothesis), position = position_dodge(0.9),vjust = -0.3, size = 2) +
  labs(x="") + scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,3650)
dgls <- ggplot(toxSupport, aes(x=Support,y=dGLStoxCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = dGLStoxCounts, group = Hypothesis), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") + scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,3650)
ggi <- ggplot(toxSupport, aes(x=Support,y=GGItoxCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGItoxCounts, group = Hypothesis), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") + scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,3650)
ggir1 <- ggplot(toxSupport, aes(x=Support,y=GGItoxSupportR1)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGItoxSupportR1, group = Hypothesis), position = position_dodge(0.9),vjust = -0.3, size = 2)+
  labs(x="") + scale_color_manual(values=myList) + scale_fill_manual(values=myList) + 
  theme_bw() + theme(panel.border = element_blank()) + ylim(0,3650)

a95 <- ggarrange(bf,dgls,ggi, ncol=1, nrow=3, common.legend = TRUE, legend="right") 
ar1 <- ggarrange(bf,dgls,ggir1, ncol=1, nrow=3, common.legend = TRUE, legend="right")

#ggsave("Singhal_SupportBarR1_tox.jpg", plot=ar1, width = 20, height = 20, units = "cm")
#ggsave("Singhal_SupportBar95_tox.jpg", plot=a95, width = 20, height = 30, units = "cm")

