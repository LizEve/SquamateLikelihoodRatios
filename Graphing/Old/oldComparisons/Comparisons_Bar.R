rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("/Users/ChatNoir/Projects/Squam/Graphs")


#load("/Users/ChatNoir/Projects/Squam/Graphs/BF/Singhal_BFdata.RData")
#load("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/Singhal_dGLSdata.RData")
#load("/Users/ChatNoir/Projects/Squam/Graphs/GGI/Singhal_GGIdata.RData")

load("/Users/ChatNoir/Projects/Squam/Graphs/BF/Streicher_BFdata.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/dGLS/Streicher_dGLSdata.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/GGI/Streicher_GGIdata.RData")

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

####### Caclulate strong support 

# BF - Count strong support loci 
BFallSupport <- allBF %>% group_by(Hypothesis,Support) %>% dplyr::summarise(BFallCounts=n())
BFallSupport$Support <- factor(BFallSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))
BFtoxSupport <- toxBF %>% group_by(Hypothesis,Support) %>% dplyr::summarise(BFtoxCounts=n())
BFtoxSupport$Support <- factor(BFtoxSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))

BFpairSupport <- pairBF %>% group_by(Hypothesis,Support) %>% dplyr::summarise(BFpairCounts=n())
BFpairSupport$Support <- factor(BFpairSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))


# dGLS - Count strong support loci 

dGLSallSupport <- alldGLS %>% group_by(Hypothesis,Support) %>% dplyr::summarise(dGLSallCounts=n()) 
dGLSallSupport$Support <- factor(dGLSallSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))

dGLStoxSupport <- toxdGLS %>% group_by(Hypothesis,Support) %>% dplyr::summarise(dGLStoxCounts=n())
dGLStoxSupport$Support <- factor(dGLStoxSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))

dGLSpairSupport <- pairdGLS %>% group_by(Hypothesis,Support) %>% dplyr::summarise(dGLSpairCounts=n())
dGLSpairSupport$Support <- factor(dGLSpairSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))

# GGI - Count strong support loci. rank 1 and 95

GGIallSupport <- allAU %>% group_by(Hypothesis,Support) %>% dplyr::summarise(GGIallCounts=n())
GGIallSupport$Support <- factor(GGIallSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))
GGIallSupportAU <- allAU %>% group_by(Hypothesis,SupportAU) %>% dplyr::summarise(GGIallAUCounts=n()) %>%
  rename(Support=SupportAU)
GGIallSupportAU$Support <- factor(GGIallSupportAU$Support, levels=c("Strong_Against","Ambiguous","Strong")) 

GGItoxSupport <- toxAU %>% group_by(Hypothesis,Support) %>% dplyr::summarise(GGItoxCounts=n())
GGItoxSupport$toxSupport <- factor(GGItoxSupport$Support, levels=c("Strong_Against","Ambiguous","Strong"))
GGItoxSupportAU <- toxAU %>% group_by(Hypothesis,SupportAU) %>% dplyr::summarise(GGItoxAUCounts=n()) %>%
  rename(Support=SupportAU)
GGItoxSupportAU$toxSupport <- factor(GGItoxSupportAU$Support, levels=c("Strong_Against","Ambiguous","Strong"))

# Pairwise GGI  
GGIpairSupportAI <- aiSAU %>%  group_by(Hypothesis,Support) %>% dplyr::summarise(GGIpairCountsAI=n())
GGIpairSupportAI$Support <- factor(GGIpairSupportAI$Support, levels=c("Ambiguous","Strong"))
GGIpairSupportAUAI <- aiSAU %>%  group_by(Hypothesis,SupportAU) %>% dplyr::summarise(GGIpairCountsAUAI=n())
GGIpairSupportAUAI$pairSupport <- factor(GGIpairSupportAUAI$SupportAU, levels=c("Strong_Against","Ambiguous","Strong"))

GGIpairSupportSA <- saSAU %>%  group_by(Hypothesis,Support) %>% dplyr::summarise(GGIpairCountsSA=n())
GGIpairSupportSA$Support <- factor(GGIpairSupportSA$Support, levels=c("Ambiguous","Strong"))
GGIpairSupportAUSA <- saSAU %>%  group_by(Hypothesis,SupportAU) %>% dplyr::summarise(GGIpairCountsAUSA=n())
GGIpairSupportAUSA$pairSupport <- factor(GGIpairSupportAUSA$SupportAU, levels=c("Strong_Against","Ambiguous","Strong"))

GGIpairSupportSI <- siSAU %>%  group_by(Hypothesis,Support) %>% dplyr::summarise(GGIpairCountsSI=n())
GGIpairSupportSI$Support <- factor(GGIpairSupportSI$Support, levels=c("Ambiguous","Strong"))
GGIpairSupportAUSI <- siSAU %>%  group_by(Hypothesis,SupportAU) %>% dplyr::summarise(GGIpairCountsAUSI=n())
GGIpairSupportAUSI$pairSupport <- factor(GGIpairSupportAUSI$SupportAU, levels=c("Strong_Against","Ambiguous","Strong"))


### Merge all into one file 

allSupport <- merge(GGIallSupportAU,
                    merge(GGIallSupport,
                          merge(BFallSupport,
                                dGLSallSupport,by=c("Hypothesis","Support")),
                          by=c("Hypothesis","Support")),
                    by=c("Hypothesis","Support"))
allSupport$Support <- factor(allSupport$Support, levels=c("Strong","Ambiguous","Strong_Against"))
allSupport <- allSupport %>% 
  mutate(Hypothesis = str_replace(Hypothesis, "ToxSnIg", "Snake Iguania")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "ToxSnAng", "Snake Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "ToxAngIg", "Iguania Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "Sclero", "Scleroglossa ")) 
allSupport$Hypothesis <- factor(allSupport$Hypothesis, 
                                levels= c("Scleroglossa ","Iguania Anguimorph", "Snake Anguimorph", "Snake Iguania"))


toxSupport <- merge(GGItoxSupportAU,
                    merge(GGItoxSupport,
                          merge(BFtoxSupport,
                                dGLStoxSupport,by=c("Hypothesis","Support")),
                          by=c("Hypothesis","Support")),
                    by=c("Hypothesis","Support"))
toxSupport$Support <- factor(toxSupport$Support, levels=c("Strong","Ambiguous","Strong_Against"))
toxSupport <- toxSupport %>% 
  mutate(Hypothesis = str_replace(Hypothesis, "ToxSnIg", "Snake Iguania")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "ToxSnAng", "Snake Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "ToxAngIg", "Iguania Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "Sclero", "Scleroglossa ")) 
toxSupport$Hypothesis <- factor(toxSupport$Hypothesis, 
                                levels= c("Scleroglossa ","Iguania Anguimorph", "Snake Anguimorph", "Snake Iguania"))
############# Bar graphs Pair BF

pairSupport <- BFpairSupport
strong <- pairSupport %>% filter(Support != "Ambiguous") %>% rename(PairTest=Hypothesis)
strong <- strong %>% ungroup() %>% group_by(PairTest)%>% mutate(Hypothesis=PairTest) %>%
  mutate(Hypothesis = case_when(Support == "Strong" ~ "Scleroglossa",
                                Support == "Strong_Against" ~ Hypothesis)) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "SvSA", "Snake Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "SvSI", "Snake Iguania")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "SvAI", "Iguania Anguimorph")) 
strong$Hypothesis <- factor(strong$Hypothesis, 
                                 levels= c("Scleroglossa","Iguania Anguimorph", "Snake Anguimorph", "Snake Iguania"))



levels(strong$Hypothesis) <- gsub(" ", "\n", levels(strong$Hypothesis))
quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
c1 <- c("orange","#38598CFF")
p1 <- strong %>% filter(PairTest == "SvAI")
p2 <- strong %>% filter(PairTest == "SvAI")
bf1 <- ggplot(p1, aes(x=Hypothesis,y=BFpairCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = BFpairCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c1) +
  labs(x="", y="BF > 10", size=24) + 
  scale_color_manual(values=c1) + 
  scale_fill_manual(values=c1) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,1800)
bf1 
ggsave("Streicher_Bar_pair_StrongSupport_BF1.pdf", plot=bf1, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

c2 <- c("orange","#2BB07FFF")
p2 <- strong %>% filter(PairTest == "SvSA")
bf2 <- ggplot(p2, aes(x=Hypothesis,y=BFpairCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = BFpairCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c2) +
  labs(x="", y="BF > 10", size=24) + 
  scale_color_manual(values=c2) + 
  scale_fill_manual(values=c2) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,1800)
bf2 
ggsave("Streicher_Bar_pair_StrongSupport_BF2.pdf", plot=bf2, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

c3 <- c("orange","#C2DF23FF")
p3 <- strong %>% filter(PairTest == "SvSI")
bf3 <- ggplot(p3, aes(x=Hypothesis,y=BFpairCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = BFpairCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c3) +
  labs(x="", y="BF > 10", size=24) + 
  scale_color_manual(values=c3) + 
  scale_fill_manual(values=c3) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,1800)
bf3 
ggsave("Streicher_Bar_pair_StrongSupport_BF3.pdf", plot=bf3, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

############# Bar graphs Pair dGLS

pairSupport <- dGLSpairSupport
strong <- pairSupport %>% filter(Support != "Ambiguous") %>% rename(PairTest=Hypothesis)
strong <- strong %>% ungroup() %>% group_by(PairTest)%>% mutate(Hypothesis=PairTest) %>%
  mutate(Hypothesis = case_when(Support == "Strong" ~ "Scleroglossa",
                                Support == "Strong_Against" ~ Hypothesis)) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "SvSA", "Snake Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "SvSI", "Snake Iguania")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "SvAI", "Iguania Anguimorph")) 
strong$Hypothesis <- factor(strong$Hypothesis, 
                            levels= c("Scleroglossa","Iguania Anguimorph", "Snake Anguimorph", "Snake Iguania"))



levels(strong$Hypothesis) <- gsub(" ", "\n", levels(strong$Hypothesis))
quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
c1 <- c("orange","#38598CFF")
p1 <- strong %>% filter(PairTest == "SvAI")
dgls1 <- ggplot(p1, aes(x=Hypothesis,y=dGLSpairCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = dGLSpairCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c1) +
  labs(x="", y="dGLS > 10", size=24) + 
  scale_color_manual(values=c1) + 
  scale_fill_manual(values=c1) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,2350)
dgls1 
ggsave("Streicher_Bar_pair_StrongSupport_dGLS1.pdf", plot=dgls1, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

c2 <- c("orange","#2BB07FFF")
p2 <- strong %>% filter(PairTest == "SvSA")
dgls2 <- ggplot(p2, aes(x=Hypothesis,y=dGLSpairCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = dGLSpairCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c2) +
  labs(x="", y="dGLS > 10", size=24) + 
  scale_color_manual(values=c2) + 
  scale_fill_manual(values=c2) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,2350)
dgls2 
ggsave("Streicher_Bar_pair_StrongSupport_dGLS2.pdf", plot=dgls2, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

c3 <- c("orange","#C2DF23FF")
p3 <- strong %>% filter(PairTest == "SvSI")
dgls3 <- ggplot(p3, aes(x=Hypothesis,y=dGLSpairCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = dGLSpairCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c3) +
  labs(x="", y="dGLS > 10", size=24) + 
  scale_color_manual(values=c3) + 
  scale_fill_manual(values=c3) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,2350)
dgls3 
ggsave("Streicher_Bar_pair_StrongSupport_dGLS3.pdf", plot=dgls3, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

############# Bar graphs Pair GGI

pairSupport <- GGIpairSupportAI
p1 <- pairSupport %>% filter(Support != "Ambiguous") %>% 
  rename(PairTest=Hypothesis) %>%
  mutate(Hypothesis=PairTest) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "ToxAngIg", "Iguania Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "Sclero", "Scleroglossa"))
p1$Hypothesis <- factor(p1$Hypothesis,levels= c("Scleroglossa","Iguania Anguimorph"))
levels(p1$Hypothesis) <- gsub(" ", "\n", levels(p1$Hypothesis))

pairSupport <- GGIpairSupportSA
p2 <- pairSupport %>% filter(Support != "Ambiguous") %>% 
  rename(PairTest=Hypothesis) %>%
  mutate(Hypothesis=PairTest) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "ToxSnAng", "Snake Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "Sclero", "Scleroglossa"))
p2$Hypothesis <- factor(p2$Hypothesis,levels= c("Scleroglossa","Snake Anguimorph"))
levels(p2$Hypothesis) <- gsub(" ", "\n", levels(p2$Hypothesis))

pairSupport <- GGIpairSupportSI
p3 <- pairSupport %>% filter(Support != "Ambiguous") %>% 
  rename(PairTest=Hypothesis) %>%
  mutate(Hypothesis=PairTest) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "ToxSnIg", "Snake Iguania")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "Sclero", "Scleroglossa"))
p3$Hypothesis <- factor(p3$Hypothesis,levels= c("Scleroglossa","Snake Iguania"))
levels(p3$Hypothesis) <- gsub(" ", "\n", levels(p3$Hypothesis))

quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")
c1 <- c("orange","#38598CFF")
c2 <- c("orange","#2BB07FFF")
c3 <- c("orange","#C2DF23FF")
ggi1 <- ggplot(p1, aes(x=Hypothesis,y=GGIpairCountsAI)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGIpairCountsAI, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c1) +
  labs(x="", y="GGI rank 1", size=24) + 
  scale_color_manual(values=c1) + 
  scale_fill_manual(values=c1) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,2600)
ggi1 
ggsave("Streicher_Bar_pair_StrongSupport_GGI1.pdf", plot=ggi1, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

ggi2 <- ggplot(p2, aes(x=Hypothesis,y=GGIpairCountsSA)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGIpairCountsSA, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c2) +
  labs(x="", y="GGI rank 1", size=24) + 
  scale_color_manual(values=c2) + 
  scale_fill_manual(values=c2) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,2600)
ggi2 
ggsave("Streicher_Bar_pair_StrongSupport_GGI2.pdf", plot=ggi2, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

ggi3 <- ggplot(p3, aes(x=Hypothesis,y=GGIpairCountsSI)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGIpairCountsSI, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=c3) +
  labs(x="", y="GGI rank 1", size=24) + 
  scale_color_manual(values=c3) + 
  scale_fill_manual(values=c3) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,2600)
ggi3
ggsave("Streicher_Bar_pair_StrongSupport_GGI3.pdf", plot=ggi3, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

############# Bar graphs ALL ############# 


strong <- allSupport[allSupport$Support == "Strong", ] 
levels(strong$Hypothesis) <- gsub(" ", "\n", levels(strong$Hypothesis))
quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")

bf <- ggplot(strong, aes(x=Hypothesis,y=BFallCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = BFallCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=myList) +
  labs(x="", y="BF > 10", size=24) + 
  scale_color_manual(values=myList) + 
  scale_fill_manual(values=myList) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
        ) + 
  ylim(0,200)
bf 
ggsave("Streicher_Bar_StrongSupport_BF.pdf", plot=bf, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

dgls <- ggplot(strong, aes(x=Hypothesis,y=dGLSallCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = dGLSallCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=myList) +
  labs(x="", y="dGLS > 0.5", size=24) + 
  scale_color_manual(values=myList) + 
  scale_fill_manual(values=myList) +
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
  ) + 
  ylim(0,2200)
ggsave("Streicher_Bar_StrongSupport_dGLS.pdf", plot=dgls,
       width = 9.5, height = 8, units = "in", device = 'pdf',
       bg = "transparent")

ggi <- ggplot(strong, aes(x=Hypothesis,y=GGIallCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGIallCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=myList) +
  labs(x="", y="GGI rank 1", size=24) + 
  scale_color_manual(values=myList) + 
  scale_fill_manual(values=myList) +
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
  ) + 
  ylim(0,1600)
ggsave("Streicher_Bar_StrongSupport_GGI.pdf", plot=ggi, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

# same pattern 
ggi95 <- ggplot(strong, aes(x=Hypothesis,y=GGIallAUCounts)) +
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGIallAUCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=myList) +
  labs(x="", y="GGI p > 0.95", size=24) + 
  scale_color_manual(values=myList) + 
  scale_fill_manual(values=myList) +
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
  ) + 
  ylim(0,1700)
ggsave("Streicher_Bar_StrongSupport_GGI95.pdf", plot=ggi95, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")


############# Bar graphs TOX ############# 


strong <- toxSupport[toxSupport$Support == "Strong", ] 
levels(strong$Hypothesis) <- gsub(" ", "\n", levels(strong$Hypothesis))
quartz()
myList <- c("#38598CFF","#2BB07FFF","#C2DF23FF")

bf <- ggplot(strong, aes(x=Hypothesis,y=BFtoxCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = BFtoxCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=myList) +
  labs(x="", y="BF > 10", size=24) + 
  scale_color_manual(values=myList) + 
  scale_fill_manual(values=myList) +
  theme_classic() + 
  theme(axis.text.x = element_text(size=24, color="black"),
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
  ) + 
  ylim(0,1700)
bf 
ggsave("Streicher_Bar_Tox_StrongSupport_BF.pdf", plot=bf, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

dgls <- ggplot(strong, aes(x=Hypothesis,y=dGLStoxCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = dGLStoxCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=myList) +
  labs(x="", y="dGLS > 0.5", size=24) + 
  scale_color_manual(values=myList) + 
  scale_fill_manual(values=myList) +
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
  ) + 
  ylim(0,1700)
ggsave("Streicher_Bar_Tox_StrongSupport_dGLS.pdf", plot=dgls,
       width = 9.5, height = 8, units = "in", device = 'pdf',
       bg = "transparent")

ggi <- ggplot(strong, aes(x=Hypothesis,y=GGItoxCounts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGItoxCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=myList) +
  labs(x="", y="GGI rank 1", size=24) + 
  scale_color_manual(values=myList) + 
  scale_fill_manual(values=myList) +
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
  ) + 
  ylim(0,1700)
ggsave("Streicher_Bar_Tox_StrongSupport_GGI.pdf", plot=ggi, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

# same pattern 
ggi95 <- ggplot(strong, aes(x=Hypothesis,y=GGItoxAUCounts)) +
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = GGItoxAUCounts, group = Hypothesis), 
            position = position_dodge(0.9),vjust = -0.3, size = 10, color=myList) +
  labs(x="", y="GGI p > 0.95", size=24) + 
  scale_color_manual(values=myList) + 
  scale_fill_manual(values=myList) +
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
  ) + 
  ylim(0,1700)
ggsave("Streicher_Bar_Tox_StrongSupport_GGI95.pdf", plot=ggi95, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")


