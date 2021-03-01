rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# set dataset 

dataset <- "Streicher"

# AI
#h1 <- "AIvSA"
#h2 <- "AIvSI"
#h3 <- "AIvAvg"
## SA
#h <- "SAdGLS"
#h1 <- "SAvSI"
#h2 <- "SAvSI"
#h3 <- "SAvAvg"
## SI
#h <- "SIdGLS"
#h1 <- "SIvAI"
#h2 <- "SIvSA"
#h3 <- "SIvAvg"



setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))


V <- function(df.g,xval.g,yval.g,ylablab.g,ytic.g,hline.g,cc.g,title.g){
  GL_violin <- ggplot(df.g, aes(x=xval.g, y=yval.g, fill=xval.g)) + 
    geom_violin(trim=TRUE) +  
    geom_point(shape = 18,size=0.5, position = position_jitterdodge(), color='black',alpha=1)+
    scale_color_manual(values=cc.g) + scale_fill_manual(values=cc.g) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=10, color="black"),
          text = element_text(size=14),
          legend.position = "none") +
    labs(y=ylablab.g,x="")  + 
    ggtitle(title.g) +
    coord_cartesian(ylim=ytic.g) +
    scale_y_continuous(breaks = ytic.g) + 
    geom_hline(yintercept=hline.g,color=c("black"), linetype="dashed", size=0.5)
  return(GL_violin)
}



# Colors 
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "#C2DF23FF"

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS"))
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))


df.1 <- mLGL
xval.1 <- mLGL$supportType
ylablab.1 <- "dGLS value"
hline.1 <- c(-0.5,0.5)

support <- "dGLS"

h1 <- "AI"
yval.1 <- mLGL$AIvSA
yval.2 <- mLGL$AIvSI
cc.1 <- c(color_AI)

h2 <- "SA"
yval.3 <- mLGL$AIvSA*(-1) # reversing values 
yval.4 <- mLGL$SAvSI
cc.2 <- c(color_SA)

h3 <- "SI"
yval.5 <- mLGL$AIvSI*(-1) # reversing values 
yval.6 <- mLGL$SAvSI*(-1) # reversing
cc.3 <- c(color_SI)


max(c(abs(yval.1),abs(yval.2),abs(yval.3),abs(yval.4),abs(yval.5),abs(yval.6)))
limit.1 <- 45
ytic.1 <- c(seq(-limit.1,limit.1,5))


title.1 <- paste(h1,h2,sep="v") 
title.2 <- paste(h1,h3,sep="v") 
title.3 <- paste(h2,h1,sep="v")  
title.4 <- paste(h2,h3,sep="v") 
title.5 <- paste(h3,h1,sep="v") 
title.6 <- paste(h3,h2,sep="v")  

GL.1 <- V(df.1,xval.1,yval.1,ylablab.1,ytic.1,hline.1,cc.1,title.1)
GL.2 <- V(df.1,xval.1,yval.2,ylablab.1,ytic.1,hline.1,cc.1,title.2)
GL.3 <- V(df.1,xval.1,yval.3,ylablab.1,ytic.1,hline.1,cc.2,title.3)
GL.4 <- V(df.1,xval.1,yval.4,ylablab.1,ytic.1,hline.1,cc.2,title.4)
GL.5 <- V(df.1,xval.1,yval.5,ylablab.1,ytic.1,hline.1,cc.3,title.5)
GL.6 <- V(df.1,xval.1,yval.6,ylablab.1,ytic.1,hline.1,cc.3,title.6)

GL <- ggarrange(GL.1, GL.2, GL.3, GL.4, GL.5, GL.6, ncol=6, nrow=1, align="h")
GL

#ggsave(paste(dataset,"_violin_scaled_",h,".pdf",sep=""), plot=GL,width = 6, height = 7, units = "in", device = 'pdf',bg = "transparent")
ggsave(paste(dataset,"_violin_wide_","all",".pdf",sep=""), plot=GL,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")

# Remove variables im gonna duplicate for BF 
rm(h, df.1,xval.1,yval.1,yval.2,yval.3,ylablab.1,ytic.1,hline.1,title.1,title.2,title.3)

df.1 <- mLBF
xval.1 <- mLBF$supportType
ylablab.1 <- "2ln(BF) value"
hline.1 <- c(-10,10)

title.1 <- paste(dataset," - ",h1,sep="") 
title.2 <- paste(dataset," - ",h2,sep="") 
title.3 <- paste(dataset," - ",h3,sep="") 

yval.1 <- mLBF$AIvSA
cc.1 <- c(color_AI)

yval.2 <- mLBF$AIvSI
cc.2 <- c(color_AI)

yval.3 <- mLBF$AIvAvg
cc.3 <- c(color_AI)

limit.1 <- 5 + round(max(abs(min(yval.1)),abs(max(yval.1))),-1)
limit.2 <- 5 + round(max(abs(min(yval.2)),abs(max(yval.2))),-1)
limit.3 <- 5 + round(max(abs(min(yval.3)),abs(max(yval.3))),-1)
max(c(limit.1,limit.2,limit.3))
limit.1 <- 100
ytic.1 <- c(seq(-limit.1,limit.1,10))

BF.1 <- V(df.1,xval.1,yval.1,ylablab.1,ytic.1,hline.1,cc.1,title.1)
BF.2 <- V(df.1,xval.1,yval.2,ylablab.1,ytic.1,hline.1,cc.2,title.2)
BF.3 <- V(df.1,xval.1,yval.3,ylablab.1,ytic.1,hline.1,cc.3,title.3)

BF <- ggarrange(BF.1, BF.2, BF.3, ncol=3, nrow=1, align="h")
BF

h <- "AIBF"
ggsave(paste(dataset,"_violin_scaled_",h,".pdf",sep=""), plot=BF,width = 6, height = 7, units = "in", device = 'pdf',bg = "transparent")


both <- ggarrange(GL.1,BF.1,GL.2,BF.2, ncol=4,nrow=1,align="h")
both
ggsave(paste(dataset,"_violin_scaled_",h3,"_BFdGLS.pdf",sep=""), plot=BF,width = 6, height = 7, units = "in", device = 'pdf',bg = "transparent")

