rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# set dataset 

dataset <- "Singhal"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Colors 
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "#C2DF23FF"


# you'll have to graph first, then reset this 
maxloci <- 1200
ymax <- seq(0,maxloci,100)

#Streicher 350,50
#Singhal 1200,100
#quartz()

# Set title and input data 
title.g <- paste(dataset," - Tox vs Sclero",sep="") 
cc.g <- color_TP
df.g <- mLGL
xval.g <- mLGL$TvS
xlablab.g <- 'dGLS values'

# Get max min for graph to set x axis values
limit.g <- 5 + round(max(abs(min(xval.g)),abs(max(xval.g))),-1)
limit.g
tic.g <- seq(-limit.g,limit.g,10)
lines.g <- c(0.5,-0.5)


GL_hist <- ggplot(data=df.g, aes(x=xval.g)) + 
  geom_histogram(binwidth = limit.g*0.01, alpha=1, position="dodge", color=cc.g, fill=cc.g)+ 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text = element_text(size=10, color="black"),
        text = element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(y="Number of Loci",x=xlablab.g)  + 
  ggtitle(title.g) +
  coord_cartesian(ylim=ymax, xlim = tic.g) +
  scale_x_continuous(breaks = c(0,tic.g)) +
  scale_y_continuous(breaks = ymax) + 
  # geom_vline(xintercept=lines.g,color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(0),color=c("black"), linetype="dashed", size=0.2)
GL_hist

  
  



# Set title and input data 
title <- ""
cc <- color_TP
df <- mLBF
xval <- mLBF$TvS
xlablab <- '2ln(BF) values'

# Get max min for graph to set x axis values
limit <- 10 + round(max(abs(min(xval)),abs(max(xval))),-1)
limit
tic <- seq(-limit,limit,20)
lines <- c(10,-10)


BF_hist <- ggplot(data=df, aes(x=xval)) + 
  geom_histogram(binwidth = limit*0.01, alpha=1, position="dodge", color=cc, fill=cc)+ 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text = element_text(size=10, color="black"),
        text = element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(y="Number of Loci",x=xlablab)  + 
  ggtitle(title) +
  coord_cartesian(ylim=ymax, xlim = tic) +
  scale_x_continuous(breaks = c(0,tic)) +
  scale_y_continuous(breaks = ymax) + 
  #geom_vline(xintercept=lines,color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(0),color=c("black"), linetype="dashed", size=0.2)


BF_hist




graph_combined <- ggarrange(GL_hist, BF_hist, ncol=1, nrow=2, align="v")
graph_combined

ggsave(paste(dataset,"_histo.pdf",sep=""), plot=graph_combined,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")


