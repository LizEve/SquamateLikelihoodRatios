rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# set dataset 

dataset1 <- "Singhal"
dataset2 <- "SinghalOG"
dataset <- "SinghalOG"

h <- "TvS"



setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data and rename data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset1,"/Calcs_",dataset1,".RData", sep=''))
mLBF1 <- mLBF
mLGL1 <- mLGL
rm(mLBF,mLGL)

load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset2,"/Calcs_",dataset2,".RData", sep=''))
mLBF2 <- mLBF
mLGL2 <- mLGL
rm(mLBF,mLGL)


## Number of overlapping
length(intersect(mLBF1$Locus,mLBF2$Locus))

# Remove any other from the larger dataset
mLBF1 <- mLBF1[mLBF1$Locus %in% intersect(mLBF1$Locus,mLBF2$Locus),]
mLBF2 <- mLBF2[mLBF2$Locus %in% intersect(mLBF2$Locus,mLBF1$Locus),]

mLGL1 <- mLGL1[mLGL1$Locus %in% intersect(mLGL1$Locus,mLGL2$Locus),]
mLGL2 <- mLGL2[mLGL2$Locus %in% intersect(mLGL2$Locus,mLGL1$Locus),]




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
maxloci <- 600
ymax <- seq(0,maxloci,100)

# Set title and input data 
title.1 <- paste(dataset1," - ",h,sep="") 
title.2 <- paste(dataset2," - ",h,sep="") 
cc.1 <- color_TP
# GL
#df.1 <- mLGL1
#df.2 <- mLGL2
#xval.1 <- mLGL1$TvS
#xval.2 <- mLGL2$TvS

# GL
df.1 <- mLBF1
df.2 <- mLBF2
xval.1 <- mLBF1$TvS
xval.2 <- mLBF2$TvS


xlablab.1 <- 'dGLS values'

xlablab.1 <- '2ln(BF) values'

# Get max min for graph to set x axis values
limit.1 <- 5 + round(max(abs(min(xval.1)),abs(max(xval.1))),-1)
limit.2 <- 5 + round(max(abs(min(xval.2)),abs(max(xval.2))),-1)
max(c(limit.1,limit.2))
limit.1 <- 210
tic.1 <- seq(-limit.1,limit.1,20)
lines.1 <- c(0.5,-0.5)

#quartz()
GL.1 <- GL(df.1,xval.1,xlablab.1,tic.1,lines.1,limit.1,cc.1,title.1)
GL.2 <- GL(df.2,xval.2,xlablab.1,tic.1,lines.1,limit.1,cc.1,title.2)

G.all <- ggarrange(GL.1,GL.2, ncol=1, nrow=2, align="v")
#quartz()
G.all

# For tox add hypoth you mentioned at the beginning 
ggsave(paste(dataset,"_histo_",h,"bf.pdf",sep=""), plot=G.all,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")

