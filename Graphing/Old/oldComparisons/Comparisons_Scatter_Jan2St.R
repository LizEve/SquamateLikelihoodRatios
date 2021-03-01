rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("/Users/ChatNoir/Projects/Squam/Graphs")

#load("/Users/ChatNoir/Projects/Squam/Graphs/Burbrink/BFcalcs_Burbrink.RData")
#load("/Users/ChatNoir/Projects/Squam/Graphs/Burbrink/dGLScalcs_Burbrink.RData")

load("/Users/ChatNoir/Projects/Squam/Graphs/Streicher/BFcalcs_Streicher.RData")
load("/Users/ChatNoir/Projects/Squam/Graphs/Streicher/dGLScalcs_Streicher.RData")
# doesnt work, no idea why 
splot <- function(df,xv,yv){
  z <- ggplot(df, aes(x=xv,y=xv)) + 
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
  return(z)
    }


######## ALL AVG Comparison #############

allBF <- allBF_St
alldGLS <- alldGLS_St
all <- merge(allBF, alldGLS, by = c("Locus","Hypothesis"))

SC <- filter(all,Hypothesis == "Sclero")
TP <- filter(all,Hypothesis == "ToxPoly")
AI <- filter(all,Hypothesis == "ToxAngIg")
SA <- filter(all,Hypothesis == "ToxSnAng")
SI <- filter(all,Hypothesis == "ToxSnIg")


quartz()
z
myList <- c("orange","#108001","#38598CFF","#2BB07FFF","#C2DF23FF")

cc <- c("orange") # Sclero

# color by locus rank 
z <- ggplot(SC, aes(x=BF,y=dGLS)) + 
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
ggsave("Streicher_Scatter_SC.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")

cc <- c("#108001") # ToxPoly

# color by locus rank 
z <- ggplot(TP, aes(x=BF,y=dGLS)) + 
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
ggsave("Streicher_Scatter_TP.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")



cc <- c("#38598CFF") # AI

# color by locus rank 
z <- ggplot(AI, aes(x=BF,y=dGLS)) + 
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
ggsave("Streicher_Scatter_AI.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")


######## TOX ONLY AVG Comparison #############

allBF <- toxBF_St
alldGLS <- toxdGLS_St
all <- merge(allBF, alldGLS, by = c("Locus","Hypothesis"))

AI <- filter(all,Hypothesis == "ToxAngIg")
SA <- filter(all,Hypothesis == "ToxSnAng")
SI <- filter(all,Hypothesis == "ToxSnIg")


cc <- c("#38598CFF") # AI

# color by locus rank 
z <- ggplot(AI, aes(x=BF,y=dGLS)) + 
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
ggsave("Streicher_Scatter_AItox.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")

######## PAIRWISE Comparison #############

allBF <- pairBF_St
alldGLS <- pairdGLS_St
all <- merge(allBF, alldGLS, by = c("Locus","Hypothesis"))


TvS <- filter(all,Hypothesis == "TvS")
AIvS <- filter(all,Hypothesis == "AIvS")
SAvS <- filter(all,Hypothesis == "SAvS")
SIvS <- filter(all,Hypothesis == "SIvS")


quartz()
z
myList <- c("orange","#108001","#38598CFF","#2BB07FFF","#C2DF23FF")

cc <- c("#108001") # TvS

# color by locus rank 
z <- ggplot(TvS, aes(x=BF,y=dGLS)) + 
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
ggsave("Streicher_Scatter_TvS.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")

cc <- c("#38598CFF") # AIvS
z <- ggplot(AIvS, aes(x=BF,y=dGLS)) + 
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
ggsave("Streicher_Scatter_AIvS.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")

cc <- c("#2BB07FFF") # SAvS
z <- ggplot(SAvS, aes(x=BF,y=dGLS)) + 
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
ggsave("Streicher_Scatter_SAvS.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")

# color by locus rank 
cc <- c("#C2DF23FF") # SIvS
z <- ggplot(SIvS, aes(x=BF,y=dGLS)) + 
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
ggsave("Streicher_Scatter_SIvS.pdf", plot=z, 
       width = 7, height = 6, units = "in", 
       device = 'pdf',bg = "transparent")

