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

all <- merge(merge(allBF, alldGLS, by = c("Locus","Hypothesis")),allGGI, by = c("Locus","Hypothesis"))

SC <- filter(all,Hypothesis == "Sclero") %>% filter(Support != "Ambiguous")
AI <- filter(all,Hypothesis == "ToxAngIg")
SA <- filter(all,Hypothesis == "ToxSnAng")
SI <- filter(all,Hypothesis == "ToxSnIg")

c1 <- c("#990000","#482173FF")
z <- ggplot(SC, aes(x=BF,y=dGLSavg)) + 
  geom_point(alpha=0.5, aes(color=factor(GGI.Rank))) + 
  scale_color_manual(values=c1) + 
  scale_fill_manual(values=c1) +
  labs(x="dGLS value",y="BF value",size=24) +
  theme_bw() + theme(panel.border = element_blank()) +
  theme(
    axis.text.x= element_text(size=24, color="black"),
    axis.text.y = element_text(size=24, color="black"),
    text = element_text(size=30),
    panel.grid.major = element_line(colour="grey88"),
    panel.grid.minor = element_line(colour="grey88"),
    legend.position = "top", 
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) + ylim(-50,60) + 
  labs(color = "GGI rank") 
z
ggsave("ToyExample_Scatter2.pdf", plot=z, 
       width = 7, height = 6, units = "in", device = 'pdf',
       bg = "transparent")

SC <- filter(all,Hypothesis == "Sclero")
z <- ggplot(SC, aes(x=BF,y=dGLSavg)) + 
  geom_point(size=5,alpha=0.5, aes(color=factor(GGI.Rank))) +
  labs(x="dGLS value",y="BF value",size=24) +
  theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  theme(
    axis.text.x= element_text(size=24, color="black"),
    axis.text.y = element_text(size=24, color="black"),
    text = element_text(size=30),
    panel.grid.major = element_line(colour="grey88"),
    panel.grid.minor = element_line(colour="grey88"),
    legend.position = "top", 
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.key=element_blank() # get rid of legend bg
    ) + ylim(-50,60) + 
  labs(color = "GGI rank") 
z
ggsave("ToyExample_Scatter2.pdf", plot=z, 
       width = 7, height = 6, units = "in", device = 'pdf',
       bg = "transparent")


SC <- filter(all,Hypothesis == "Sclero")
z <- ggplot(SC, aes(x=BF,y=dGLSavg)) + 
  geom_point(size=5,alpha=0.5, color="#990000") +
  labs(x="dGLS value",y="BF value",size=24) +
  theme_bw() + theme(panel.border = element_blank()) +
  geom_vline(xintercept=c(10),color=c("black"), linetype="dashed", size=0.5) +
  geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.5) +
  theme(
    axis.text.x= element_text(size=24, color="black"),
    axis.text.y = element_text(size=24, color="black"),
    text = element_text(size=30),
    panel.grid.major = element_line(colour="grey88"),
    panel.grid.minor = element_line(colour="grey88"),
    legend.position = "top", 
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.key=element_blank() # get rid of legend bg
  ) + ylim(-50,60) + 
  labs(color = "GGI rank") 
z
ggsave("ToyExample_Scatter3.pdf", plot=z, 
       width = 7, height = 6, units = "in", device = 'pdf',
       bg = "transparent")
