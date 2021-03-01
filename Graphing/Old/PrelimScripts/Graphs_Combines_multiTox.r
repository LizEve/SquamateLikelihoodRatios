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
color_SI <- "yellow4"

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS"))
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))

# Want to stack datasets 

reformatdf <- function(df,h,h1,s1,keepCols){
  # keepCols is a vector of the locus and support type column indices
  # Grab column number  
  c1 <- match(h1,names(df))
  # Grab and rename column, add hypothesis column
  a.df <- df[,append(keepCols,c1)] 
  names(a.df)[names(a.df) == h1] <- h
  a.df$Hypothesis <- rep(h1,length(a.df[,1]))
  c <- match(h,names(a.df))
  # Adjust direction of support if needed. ie if column is AIvSA, but you want to know SAvAI
  a.df[,c] <- a.df[,c]*(s1)
  return(a.df)
}

# reformat 
AI.g <- bind_rows(reformatdf(mLGL,"AI","AIvSA",1,c(1,14)),
                  reformatdf(mLGL,"AI","AIvSI",1,c(1,14)))

SA.g <- bind_rows(reformatdf(mLGL,"SA","AIvSA",-1,c(1,14)),
                  reformatdf(mLGL,"SA","SAvSI",1,c(1,14)))

SI.g <- bind_rows(reformatdf(mLGL,"SI","AIvSI",-1,c(1,14)),
                  reformatdf(mLGL,"SI","SAvSI",-1,c(1,14)))

AI.b <- bind_rows(reformatdf(mLBF,"AI","AIvSA",1,c(1,14)),
                  reformatdf(mLBF,"AI","AIvSI",1,c(1,14)))

SA.b <- bind_rows(reformatdf(mLBF,"SA","AIvSA",-1,c(1,14)),
                  reformatdf(mLBF,"SA","SAvSI",1,c(1,14)))

SI.b <- bind_rows(reformatdf(mLBF,"SI","AIvSI",-1,c(1,14)),
                  reformatdf(mLBF,"SI","SAvSI",-1,c(1,14)))

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Histogram -----------------------------------------Histogram---------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

H <- function(df.g,xval.g,xlablab.g,xtic.g,lines.g,limit.g,bin.g,cc.g,title.g,ytic.g){
  GL_hist <- ggplot(data=df.g, aes(x=xval.g)) + 
    #geom_histogram(binwidth = limit.g*0.001, alpha=1, position="dodge", color=cc.g, fill=cc.g)+
    geom_histogram(binwidth = bin.g, alpha=1, position="dodge", color=cc.g, fill=cc.g)+
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=10, color="black"),
          text = element_text(size=14),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)) +
    labs(y="Number of Loci",x=xlablab.g)  + 
    ggtitle(title.g) +
    coord_cartesian(ylim=ytic.g, xlim = xtic.g) +
    scale_x_continuous(breaks = c(0,xtic.g)) +
    scale_y_continuous(breaks = ytic.g) + 
    geom_vline(xintercept=lines.g,color=c("black"), linetype="dashed", size=0.5)
    #geom_vline(xintercept=c(0),color=c("black"), linetype="dashed", size=0.2)
  return(GL_hist)
}

# These will stay the same across all graphs 

# x lable and where lines are 
xlablab.1 <- "dGLS value"
lines.1 <- c(-0.5,0.5)

# Set max X value, and max Y value
max(c(abs(AI.g$AI),abs(SA.g$SA),abs(SI.g$SI)))
limit.1 <- 45
xtic.1 <- c(seq(-limit.1,limit.1,5))
maxloci.1 <- 100
ytic.1 <- seq(0,maxloci.1,10)
bin.1 <- 0.05 # 0.5/10 = 0.05

quartz()
AI.gg <- H(AI.g,AI.g$AI,xlablab.1,xtic.1,lines.1,limit.1,bin.1,color_AI,"AI",ytic.1)
AI.gg
SA.gg <- H(SA.g,SA.g$SA,xlablab.1,xtic.1,lines.1,limit.1,bin.1,color_SA,"SA",ytic.1)
SI.gg <- H(SI.g,SI.g$SI,xlablab.1,xtic.1,lines.1,limit.1,bin.1,color_SI,"SI",ytic.1)

GL <- ggarrange(AI.gg, SA.gg, SI.gg, ncol=1, nrow=3, align="v")
GL

#ggsave(paste(dataset,"_histo_tox_gl.pdf",sep=""), plot=GL,width = 9, height = 9, units = "in", device = 'pdf',bg = "transparent")

# These will stay the same across all graphs 

# x lable and where lines are 
xlablab.2 <- "2ln(BF) value"
lines.2 <- c(-10,10)

# Set max X value, and max Y value
max(c(abs(AI.b$AI),abs(SA.b$SA),abs(SI.b$SI)))
limit.2 <- 115
xtic.2 <- c(seq(-limit.2,limit.2,10))
maxloci.2 <- 500
ytic.2 <- seq(0,maxloci.2,100)
bin.2 <- 1 # 10/10 = 1

#quartz()
AI.bg <- H(AI.b,AI.b$AI,xlablab.2,xtic.2,lines.2,limit.2,bin.2,color_AI,"AI",ytic.2)
SA.bg <- H(SA.b,SA.b$SA,xlablab.2,xtic.2,lines.2,limit.2,bin.2,color_SA,"SA",ytic.2)
SI.bg <- H(SI.b,SI.b$SI,xlablab.2,xtic.2,lines.2,limit.2,bin.2,color_SI,"SI",ytic.2)

BF <- ggarrange(AI.bg, SA.bg, SI.bg, ncol=1, nrow=3, align="v")
BF

# ggsave(paste(dataset,"_histo_tox_bf.pdf",sep=""), plot=BF,width = 7, height = 9, units = "in", device = 'pdf',bg = "transparent")

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Violin ---------------------------------------------- Violin --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------


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


# These will stay the same across all graphs 

# x lable and where lines are 
ylablab.1 <- "dGLS value"
lines.1 <- c(-0.5,0.5)

# Set max X value, and max Y value
max(c(abs(AI.g$AI),abs(SA.g$SA),abs(SI.g$SI)))
limit.1 <- 10
ytic.1 <- c(seq(-limit.1,limit.1,1))

#quartz()
AI.gg <- V(AI.g,AI.g$supportType,AI.g$AI,ylablab.1,ytic.1,lines.1,color_AI,"AI")
SA.gg <- V(SA.g,SA.g$supportType,SA.g$SA,ylablab.1,ytic.1,lines.1,color_SA,"SA")
SI.gg <- V(SI.g,SI.g$supportType,SI.g$SI,ylablab.1,ytic.1,lines.1,color_SI,"SI")
GL <- ggarrange(AI.gg, SA.gg, SI.gg, ncol=3, nrow=1, align="h")
GL

#ggsave(paste(dataset,"_violin_scaled_tox_gl.pdf",sep=""), plot=GL,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")


# These will stay the same across all graphs 

# x lable and where lines are 
ylablab.2 <- "2ln(BF) value"
lines.2 <- c(-10,10)

# Set max X value, and max Y value
max(c(abs(AI.g$AI),abs(SA.g$SA),abs(SI.g$SI)))
limit.2 <- 40
ytic.2 <- c(seq(-limit.2,limit.2,5))

#quartz()
AI.bg <- V(AI.b,AI.b$supportType,AI.b$AI,ylablab.2,ytic.2,lines.2,color_AI,"AI")
SA.bg <- V(SA.b,SA.b$supportType,SA.b$SA,ylablab.2,ytic.2,lines.2,color_SA,"SA")
SI.bg <- V(SI.b,SI.b$supportType,SI.b$SI,ylablab.2,ytic.2,lines.2,color_SI,"SI")
BF <- ggarrange(AI.bg, SA.bg, SI.bg, ncol=3, nrow=1, align="h")
BF

#ggsave(paste(dataset,"_violin_zoom_tox_bf.pdf",sep=""), plot=BF,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")

#------------------------------------------------------------------------------------------------------------------------------------------------------

# BoxPlot ---------------------------------------------- BoxPlot --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
quartz()
#FAIL 
max(c(abs(AI.g$AI),abs(SA.g$SA),abs(SI.g$SI)))
limit.2 <- 20
ytic.2 <- c(seq(-limit.2,limit.2,1))
p <- ggplot(AI.b, aes(x=Hypothesis, y=AI)) + 
  geom_boxplot() + coord_cartesian(ylim=ytic.2) +
  scale_y_continuous(breaks = ytic.2)
p

limit.2 <- 2
ytic.2 <- c(seq(-limit.2,limit.2,0.1))
q <- ggplot(AI.g, aes(x=Hypothesis, y=AI)) + 
  geom_boxplot() + coord_cartesian(ylim=ytic.2) +
  scale_y_continuous(breaks = ytic.2)
quartz()
q



#------------------------------------------------------------------------------------------------------------------------------------------------------

# Deciles ---------------------------------------------- Deciles --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Set title and input data 
df.g <- AI.g
xval.g <- AI.g$AI
xlablab.g <- paste("dGLS - AI - ",dataset,sep="") 


df.g <- SA.g
xval.g <- SA.g$SA
xlablab.g <- paste("dGLS - SA - ",dataset,sep="") 
#quartz()


df.g <- SI.g
xval.g <- SI.g$SI
xlablab.g <- paste("dGLS - SI - ",dataset,sep="") 
#quartz()

########### BFFFF 

df.g <- AI.b
xval.g <- AI.b$AI
xlablab.g <- paste("2ln(BF) - AI - ",dataset,sep="") 


df.g <- SA.b
xval.g <- SA.b$SA
xlablab.g <- paste("2ln(BF) - SA - ",dataset,sep="") 
#quartz()


df.g <- SI.b
xval.g <- SI.b$SI
xlablab.g <- paste("2ln(BF) - SI - ",dataset,sep="") 
#quartz()

datSum.g <- ddply(df.g,.(supportType),summarise,
                  qa = quantile(xval.g,0.01),
                  qb = quantile(xval.g,0.05),
                  q1 = quantile(xval.g,0.1),
                  q2 = quantile(xval.g,0.2),
                  q3 = quantile(xval.g,0.3),
                  q4 = quantile(xval.g,0.4),
                  q5a = quantile(xval.g,0.49),
                  q5 = quantile(xval.g,0.5),
                  q5b = quantile(xval.g,0.51),
                  q6 = quantile(xval.g,0.6),
                  q7 = quantile(xval.g,0.7),
                  q8 = quantile(xval.g,0.8),
                  q9 = quantile(xval.g,0.9),
                  qy = quantile(xval.g,0.95),
                  qz = quantile(xval.g,0.99))
datSum.g$NameAlt <- 1:1
datSum.g
tics.g <- seq(-1,1,0.1) # Streicher 
tics.g <- seq(-6,6,0.5) # Streicher 
tics.g <- seq(-17,17,1) # Streicher 



GL3 <- ggplot(datSum.g,aes(ymin = NameAlt - 0.2,ymax = NameAlt + 0.2)) + 
  geom_rect(aes(xmin = qa,xmax = qz),fill = "white",colour = "black") + 
  geom_rect(aes(xmin = qb,xmax = qy),fill = "lightblue",colour = "black") + 
  geom_rect(aes(xmin = q1,xmax = q9),fill = "dodgerblue",colour = "black") + 
  geom_rect(aes(xmin = q2,xmax = q8),fill = "blue",colour = "black") + 
  geom_rect(aes(xmin = q3,xmax = q7),fill = 'midnightblue',colour = "black") +
  #geom_rect(aes(xmin = q5,xmax = q5),fill = 'midnightblue',colour = "black") +
  theme_classic() + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size=16),
        axis.text = element_text(size=10, color="black"),
        text = element_text(size=14)) + 
  labs(x=xlablab.g,y="") +
  coord_cartesian(xlim=tics.g) +
  scale_x_continuous(breaks = tics.g)
#quartz()
GL3


G <- ggarrange(GL1, GL2, GL3, ncol=1, nrow=3, align="v")
#quartz()
G

ggsave(paste(dataset,"_decile_tox_bf.pdf",sep=""), plot=BF,width = 9, height = 7, units = "in", device = 'pdf',bg = "transparent")



# https://stackoverflow.com/questions/46628958/making-a-specific-quantile-plot-in-r
# cool graph but my data is too concentrated in the middle 