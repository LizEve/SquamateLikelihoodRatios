rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# set dataset 

dataset <- "Streicher"
## AI
#h <- "AIdGLS"
#h1 <- "AIvSA"
#h2 <- "AIvSI"
#h3 <- "AIvAvg"
# SA
h <- "SAdGLS"
h1 <- "SAvSI"
h2 <- "SAvSI"
h3 <- "SAvAvg"
## SI
#h <- "SIdGLS"
#h1 <- "SIvAI"
#h2 <- "SIvSA"
#h3 <- "SIvAvg"



setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))


GL <- function(df.g,xval.g,xlablab.g,tic.g,lines.g,limit.g,cc.g,title.g){
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
  return(GL_hist)
}

# Colors 
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "#C2DF23FF"


# you'll have to graph first, then reset this 
maxloci <- 300
ymax <- seq(0,maxloci,50)

#Streicher 350,50
#Singhal 1200,100
#quartz()

# Set title and input data 
title.1 <- paste(dataset," - ",h1,sep="") 
title.2 <- paste(dataset," - ",h2,sep="") 
title.3 <- paste(dataset," - ",h3,sep="") 
cc.1 <- color_SI
df.1 <- mLGL
## AI
#xval.1 <- mLGL$AIvSA
#xval.2 <- mLGL$AIvSI
#xval.3 <- mLGL$AIvAvg
## SA
#xval.1 <- mLGL$AIvSA*(-1) # reversing values 
#xval.2 <- mLGL$SAvSI
#xval.3 <- mLGL$SAvAvg
## SI
#xval.1 <- mLGL$AIvSI*(-1) # reversing values 
#xval.2 <- mLGL$SAvSI*(-1) # reversing
#xval.3 <- mLGL$SIvAvg


xlablab.1 <- 'dGLS values'

# Get max min for graph to set x axis values
limit.1 <- 5 + round(max(abs(min(xval.1)),abs(max(xval.1))),-1)
limit.2 <- 5 + round(max(abs(min(xval.2)),abs(max(xval.2))),-1)
limit.3 <- 5 + round(max(abs(min(xval.3)),abs(max(xval.3))),-1)
max(c(limit.1,limit.2,limit.3))
limit.1 <- 25
tic.1 <- seq(-limit.1,limit.1,5)
lines.1 <- c(0.5,-0.5)

#quartz()
GL.1 <- GL(df.1,xval.1,xlablab.1,tic.1,lines.1,limit.1,cc.1,title.1)
GL.2 <- GL(df.1,xval.2,xlablab.1,tic.1,lines.1,limit.1,cc.1,title.2)
GL.3 <- GL(df.1,xval.3,xlablab.1,tic.1,lines.1,limit.1,cc.1,title.3)

G.all <- ggarrange(GL.1,GL.2,GL.3, ncol=1, nrow=3, align="v")
quartz()
G.all

#ggsave(paste(dataset,"_histo.pdf",sep=""), plot=graph_combined,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")
# For tox add hypoth you mentioned at the beginning 
ggsave(paste(dataset,"_histo_",h,".pdf",sep=""), plot=G.all,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")

# Remove variables im gonna duplicate for BF 
rm(h, df.1,xval.1,xval.2,xval.3,xlablab.1,limit.1,limit.2,limit.3,tic.1,lines.1)


# Set title and input data 
df.1 <- mLBF
## AI
#h <- "AIBF"
#xval.1 <- mLGL$AIvSA
#xval.2 <- mLGL$AIvSI
#xval.3 <- mLGL$AIvAvg
## SA
#h <- "SABF"
#xval.1 <- mLGL$AIvSA*(-1) # reversing values 
#xval.2 <- mLGL$SAvSI
#xval.3 <- mLGL$SAvAvg
## SI
#h <- "SIBF"
#xval.1 <- mLGL$AIvSI*(-1) # reversing values 
#xval.2 <- mLGL$SAvSI*(-1) # reversing
#xval.3 <- mLGL$SIvAvg



xlablab.1 <- '2ln(BF) values'

# Get max min for graph to set x axis values
limit.1 <- 5 + round(max(abs(min(xval.1)),abs(max(xval.1))),-1)
limit.2 <- 5 + round(max(abs(min(xval.2)),abs(max(xval.2))),-1)
limit.3 <- 5 + round(max(abs(min(xval.3)),abs(max(xval.3))),-1)
max(c(limit.1,limit.2,limit.3))
limit.1 <- 25
tic.1 <- seq(-limit.1,limit.1,5)
lines.1 <- c(10,10)

#quartz()
BF.1 <- GL(df.1,xval.1,xlablab.1,tic.1,lines.1,limit.1,cc.1,title.1)
BF.2 <- GL(df.1,xval.2,xlablab.1,tic.1,lines.1,limit.1,cc.1,title.2)
BF.3 <- GL(df.1,xval.3,xlablab.1,tic.1,lines.1,limit.1,cc.1,title.3)

B.all <- ggarrange(BF.1,BF.2,BF.3, ncol=1, nrow=3, align="v")
quartz()
B.all

#ggsave(paste(dataset,"_histo.pdf",sep=""), plot=graph_combined,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")
# For tox add hypoth you mentioned at the beginning 
ggsave(paste(dataset,"_histo_",h,".pdf",sep=""), plot=G.all,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")
