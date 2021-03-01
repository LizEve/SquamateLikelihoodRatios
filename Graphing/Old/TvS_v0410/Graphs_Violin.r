rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# set dataset 

dataset <- "Streicher"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Colors 
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "#C2DF23FF"

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS"))
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))


df.g <- mLGL
ylablab.g <- "dGLS value"
title.g <- ""
max(abs(min(df.g$TvS)),abs(max(df.g$TvS)))
#limit.g <- 50
#ytic.g <- c(seq(-limit.g,limit.g,5))
limit.g <- 10
ytic.g <- c(seq(-limit.g,limit.g,1))
cc <- c(color_TP)

#quartz()
GL_violin <- ggplot(df.g, aes(x=supportType, y=TvS, fill=supportType)) + 
  geom_violin(trim=TRUE) +  
  geom_point(shape = 18,size=0.75, position = position_jitterdodge(), color='black',alpha=1)+
  scale_color_manual(values=cc) + scale_fill_manual(values=cc) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text = element_text(size=10, color="black"),
        text = element_text(size=14),
        legend.position = "none") +
  labs(y=ylablab.g,x="")  + 
  ggtitle(title.g) +
  coord_cartesian(ylim=ytic.g) +
  scale_y_continuous(breaks = ytic.g) + 
  geom_hline(yintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)
  
GL_violin


df <- mLBF
ylablab <- "2ln(BF) value"
title <- ""
max(abs(min(df$TvS)),abs(max(df$TvS)))
#limit <- 1000
#ytic <- c(seq(-limit,limit,100))
limit <- 200
ytic <- c(seq(-limit,limit,20))
cc <- c(color_TP)

#quartz()
BF_violin <- ggplot(df, aes(x=supportType, y=TvS, fill=supportType)) + 
  geom_violin(trim=FALSE) +  
  geom_point(shape = 18,size=0.75, position = position_jitterdodge(), color='black',alpha=1)+
  scale_color_manual(values=cc) + scale_fill_manual(values=cc) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text = element_text(size=10, color="black"),
        text = element_text(size=14),
        legend.position = "none") +
  labs(y=ylablab,x="")  + 
  ggtitle(title) +
  coord_cartesian(ylim=ytic) +
  scale_y_continuous(breaks = ytic) + 
  geom_hline(yintercept=c(-10,10),color=c("black"), linetype="dashed", size=0.5)

BF_violin

  

graph_combined <- ggarrange(GL_violin, BF_violin, ncol=2, nrow=1, align="h")
graph_combined

ggsave(paste(dataset,"_violin","_scaled",".pdf",sep=""), plot=graph_combined,width = 4, height = 6, units = "in", device = 'pdf',bg = "transparent")

##Streicher 
#limit.g <- 50
#ytic.g <- c(seq(-limit.g,limit.g,5))
#limit.g <- 10
#ytic.g <- c(seq(-limit.g,limit.g,1))
#limit <- 1000
#ytic <- c(seq(-limit,limit,100))
#limit <- 200
#ytic <- c(seq(-limit,limit,20))

##Singhal
#limit <- 1600
#limit <- 1000
#ytic <- c(seq(-limit,limit,100))
#limit.g <-80
#limit.g <- 50
#ytic.g <- c(seq(-limit.g,limit.g,5))

##Reeder
#limit.g <- 60
#limit.g <- 40
#ytic.g <- c(seq(-limit.g,limit.g,5))
#limit <- 1200
#limit <- 800
#ytic <- c(seq(-limit,limit,100))
