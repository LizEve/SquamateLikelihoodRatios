rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("/Users/ChatNoir/Projects/Squam/Graphs")

#load("/Users/ChatNoir/Projects/Squam/Graphs/BF/Singhal_BFdata.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/BF/Streicher_BFdata.RData")
toxBF <- toxBF %>% select(-LocusNumber,-locusRank,-Percentile)%>%
  group_by(Hypothesis) %>% 
  mutate(Percentile = ntile(desc(BF),100)) %>% 
  mutate(Decile = ntile(desc(BF),10)) %>% 
  mutate(locusRank=dense_rank(desc(BF))) %>% 
  rename(BF.Percentile=Percentile,BF.Decile=Decile, 
         BF.locusRank = locusRank, BF.Support = Support)
# This will give the same rank number if the value is identicle in two cells 
#load("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/dGLSdata.avg.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/Streicher_dGLSdata.RData")
toxdGLS <- toxdGLS %>%  select(-dGLSmax,-locusRank,-Percentile) %>%
  group_by(Hypothesis) %>% 
  mutate(Percentile = ntile(desc(dGLSavg),100)) %>% 
  mutate(Decile = ntile(desc(dGLSavg),10)) %>% 
  mutate(locusRank=dense_rank(desc(dGLSavg))) %>% 
  rename(dGLS.Percentile=Percentile, dGLS.Decile=Decile, 
         dGLS.locusRank = locusRank, dGLS.Support=Supportavg)
#load("/Users/ChatNoir/Projects/Squam/Graphs/GGI/GGIdata.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/GGI/Streicher_GGIdata.RData")
toxGGI <- toxAU0.05 %>% rename(Locus=LOCUS) %>% select(-item,-gt0.95,-locusRank) %>%
  group_by(Hypothesis) %>% 
  mutate(Percentile = ntile(desc(au),100)) %>% 
  mutate(Decile = ntile(desc(au),10)) %>% 
  mutate(locusRank=dense_rank(desc(au))) %>% 
  rename(GGI.Percentile=Percentile, GGI.Decile=Decile, 
         GGI.locusRank = locusRank, GGI.Rank=rank, GGI.Support=Support)


######## Pairwise method locusRank comparison #############
#dGLS vs BF
#BF vs GGI
#dGLS vs GGI 


tox <- merge(merge(toxBF, toxdGLS, by = c("Locus","Hypothesis")),toxGGI, by = c("Locus","Hypothesis")) %>% 
  separate(Locus, c("LocusType","LocusNumber"), "-", remove = FALSE)

AI <- filter(tox,Hypothesis == "ToxAngIg")
SA <- filter(tox,Hypothesis == "ToxSnAng")
SI <- filter(tox,Hypothesis == "ToxSnIg")

quartz()
myList <- c("#38598CFF","#2BB07FFF","#C2DF23FF")

quartz()
# color by ggi rank 
a <- ggplot(tox, aes(x=BF,y=dGLSavg)) + 
  geom_point(alpha=0.5,aes(color=factor(GGI.Rank))) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5)


# color by hypothesis 
a <- ggplot(tox, aes(x=BF,y=dGLSavg)) + 
  geom_point(alpha=0.5,aes(color=tox$Hypothesis)) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  scale_color_manual(values=myList)
b <- ggplot(tox, aes(x=BF,y=au)) + 
  geom_point(alpha=0.5, aes(color=tox$Hypothesis)) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.95),color=c("black"), linetype="dashed", size=0.5) +
  scale_color_manual(values=myList)
c <- ggplot(tox, aes(x=au,y=dGLSavg)) + 
  geom_point(alpha=0.5, aes(color=tox$Hypothesis)) + theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(0.95),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  scale_color_manual(values=myList)
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="bottom", align="v",common.legend = TRUE)
x
ggsave("Streicher_Scatter_toxbyHypoth.jpg", plot=a, width = 30, height = 10, units = "cm")





a <- ggplot(tox, aes(x=BF.locusRank,y=dGLS.locusRank,color=LocusType)) + 
  geom_point(alpha=0.5) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.3) 
b <- ggplot(tox, aes(x=BF.locusRank,y=GGI.locusRank,color=LocusType)) + 
  geom_point(alpha=0.5) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.3) 
c <- ggplot(tox, aes(x=GGI.locusRank,y=dGLS.locusRank,color=LocusType)) + 
  geom_point(alpha=0.5) +
  scale_y_reverse()+scale_x_reverse()+ theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.3) 
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="bottom", align="v")
x
#ggsave("Singhal_LocusRank_tox_LocusType.jpg", plot=x, width = 30, height = 10, units = "cm")



cc <- "grey44"
a <- ggplot(tox, aes(x=BF.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
b <- ggplot(tox, aes(x=BF.locusRank,y=GGI.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse() + theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
c <- ggplot(tox, aes(x=GGI.locusRank,y=dGLS.locusRank)) + 
  geom_point(alpha=0.5,color=cc) +
  scale_y_reverse()+scale_x_reverse()+ theme_bw() + theme(panel.border = element_blank()) +
  geom_smooth(method="lm", se=FALSE, color="grey29", alpha=0.5) 
x <- ggarrange(a,b,c, ncol=3, nrow=1, legend="none", align="v")
x
ggsave("Streicher_LocusRank_tox.jpg", plot=x, width = 30, height = 10, units = "cm")
#ggsave("Singhal_LocusRank_tox.jpg", plot=x, width = 30, height = 10, units = "cm")

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
ggsave("Streicher_LocusRank_tox_AI.jpg", plot=x, width = 30, height = 10, units = "cm")
#ggsave("Singhal_LocusRank_tox_AI.jpg", plot=x, width = 30, height = 10, units = "cm")

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

ggsave("Streicher_LocusRank_tox_SA.jpg", plot=x, width = 30, height = 10, units = "cm")
#ggsave("Singhal_LocusRank_tox_SA.jpg", plot=x, width = 30, height = 10, units = "cm")


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

ggsave("Streicher_LocusRank_tox_SI.jpg", plot=x, width = 30, height = 10, units = "cm")
#ggsave("Singhal_LocusRank_tox_SI.jpg", plot=x, width = 30, height = 10, units = "cm")





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

