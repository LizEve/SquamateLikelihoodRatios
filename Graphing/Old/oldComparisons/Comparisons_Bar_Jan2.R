rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("/Users/ChatNoir/Projects/Squam/Graphs")

load("/Users/ChatNoir/Projects/Squam/Graphs/Burbrink/BFcalcs_Burbrink.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/Burbrink/dGLScalcs_Burbrink.RData")

############# Bar graphs Pair BF

pairSupport <- pairSupportBF_B
strong <- pairSupport %>% filter(Support != "Ambiguous") %>% rename(PairTest=Hypothesis)
strong <- strong %>% ungroup() %>% group_by(PairTest)%>% mutate(Hypothesis=PairTest) %>%
  mutate(Hypothesis = case_when(Support == "Strong" ~ Hypothesis,
                                Support == "Strong_Against" ~ "Scleroglossa")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "SAvS", "Snake Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "SIvS", "Snake Iguania")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "AIvS", "Iguania Anguimorph")) %>%
  mutate(Hypothesis = str_replace(Hypothesis, "TvS", "Toxicofera")) 
strong$Hypothesis <- factor(strong$Hypothesis, 
                                 levels= c("Scleroglossa","Toxicofera","Iguania Anguimorph", "Snake Anguimorph", "Snake Iguania"))



levels(strong$Hypothesis) <- gsub(" ", "\n", levels(strong$Hypothesis))
quartz()
myList <- c("orange","#38598CFF","#2BB07FFF","#C2DF23FF")

c1 <- c("orange","#108001")
p1 <- strong %>% filter(PairTest == "TvS")
bf1 <- ggplot(p1, aes(x=Hypothesis,y=counts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = counts, group = Hypothesis), 
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
  ylim(0,180)
bf1 
ggsave("Burbrink_Bar_TvSpair_BF.pdf", plot=bf1, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

c1 <- c("orange","#38598CFF")
p1 <- strong %>% filter(PairTest == "AIvS")
p2 <- strong %>% filter(Hypothesis == "Iguania Anguimorph")
bf1 <- ggplot(p1, aes(x=Hypothesis,y=counts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = counts, group = Hypothesis), 
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
  ylim(0,180)
bf1 
ggsave("Burbrink_Bar_AIvSpair_BF.pdf", plot=bf1, 
       width = 9.5, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

c2 <- c("orange","#2BB07FFF")
p2 <- strong %>% filter(PairTest == "SAvS")
bf2 <- ggplot(p2, aes(x=Hypothesis,y=counts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = counts, group = Hypothesis), 
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
  ylim(0,180)
bf2 
ggsave("Burbrink_Bar_SAvSpair_BF.pdf", plot=bf2, 
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

############# Bar graphs ALL ############# 


strong <- allSupportBF_B[allSupportBF_B$Support == "Strong", ] 
levels(strong$Hypothesis) <- gsub(" ", "\n", levels(strong$Hypothesis))
quartz()
myList <- c("orange","#108001","#38598CFF","#2BB07FFF","#C2DF23FF")

bf <- ggplot(strong, aes(x=Hypothesis,y=counts)) + 
  geom_bar(aes(color = Hypothesis, fill = Hypothesis), 
           stat = "identity", position = position_dodge()) + 
  geom_text(aes(label = counts, group = Hypothesis), 
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


