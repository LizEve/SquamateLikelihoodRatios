rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library("viridis")
# set dataset 

#dataset <- "SinghalOG"
dataset <- "Streicher"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Colors 
# https://www.w3schools.com/colors/colors_picker.asp
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS")) %>%
  mutate_if(is.numeric,round,2)
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF")) %>%
  mutate_if(is.numeric,round,2)

# Get total number of loci 
loci <- length(mLGL$Locus)

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10 - not including average calcs
names(mLBF[7:10])
mLBF[7:10] <- mLBF[7:10]/2
# for singhal OG 
#names(mLBF[4])
#mLBF[4] <- mLBF[4]/2


# Want to stack datasets so GL and BF columns have both comparisons for each hypothesis. 

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
AI.g <- bind_rows(reformatdf(mLGL,"AI","AIvSA",1,c(1,13)),
                  reformatdf(mLGL,"AI","AIvSI",1,c(1,13)))

SA.g <- bind_rows(reformatdf(mLGL,"SA","AIvSA",-1,c(1,13)),
                  reformatdf(mLGL,"SA","SAvSI",1,c(1,13)))

SI.g <- bind_rows(reformatdf(mLGL,"SI","AIvSI",-1,c(1,13)),
                  reformatdf(mLGL,"SI","SAvSI",-1,c(1,13)))

AI.b <- bind_rows(reformatdf(mLBF,"AI","AIvSA",1,c(1,13)),
                  reformatdf(mLBF,"AI","AIvSI",1,c(1,13)))

SA.b <- bind_rows(reformatdf(mLBF,"SA","AIvSA",-1,c(1,13)),
                  reformatdf(mLBF,"SA","SAvSI",1,c(1,13)))

SI.b <- bind_rows(reformatdf(mLBF,"SI","AIvSI",-1,c(1,13)),
                  reformatdf(mLBF,"SI","SAvSI",-1,c(1,13)))

TS.g <- reformatdf(mLGL,"TS","TvS",1,c(1,13)) %>% select(-SIvAvg)
TS.b <- reformatdf(mLBF,"TS","TvS",1,c(1,13)) %>% select(-SIvAvg)

# Singhal OG
#TS.g <- reformatdf(mLGL,"TS","TvS",1,c(1,3))
#TS.b <- reformatdf(mLBF,"TS","TvS",1,c(1,3))

# MERGE 

AI.a <- left_join(AI.g, AI.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
SA.a <- left_join(SA.g, SA.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
SI.a <- left_join(SI.g, SI.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
TS.a <- left_join(TS.g, TS.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))

rm(AI.g,AI.b,SA.b,SA.g,SI.b,SI.g,TS.b,TS.g)

library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_TvSLoci.xlsx",sep='')

#write.xlsx2(as.data.frame(TS.a), file=fname, sheetName="loci", row.names=FALSE)




#------------------------------------------------------------------------------------------------------------------------------------------------------

# Scatter ---------------------------------------------- Scatter BF vs dGLS --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------


## Scatter graph, colored by hypothesis 

S <- function(df,xcol,ycol,x.tic,y.tic,cc){
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]])) + 
    geom_point(alpha=1, aes(color=df$Hypothesis), size=0.5) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=12, color="black"),
      text = element_text(size=14),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x='dGLS',y='ln(BF)',color="Hypothesis") + 
    scale_color_manual(values=cc) + 
    guides(colour = guide_legend(override.aes = list(size=2))) 
  
  return(scat)
}

df <- TS.a

# Get max min for graph x = dgls y=bf
max(abs(df[[2]]))
x.lims <- 60 # 50, 45
x.ts <- seq(-x.lims,x.lims,20)
max(abs(df[[4]]))
y.lims <- 60 #1000, 120
y.ts <- seq(-y.lims,y.lims,20)

cc <- c("springgreen4")
TS.s <- S(TS.a,2,4,x.ts,y.ts,cc) + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)
quartz()
x <- TS.s + geom_point(size=1.5, color=color_TP) + theme(legend.position = "none")
x
ggsave(paste(dataset,"_scatterTS.pdf",sep=""), plot=x,width = 4, height = 3, units = "in", device = 'pdf',bg = "transparent")

df <- AI.a

# Get max min for graph x = dgls y=bf
max(abs(df[[3]]))
x.lim <- 180 # 50, 45
x.t <- seq(-x.lim,x.lim,50)
max(abs(df[[6]]))
y.lim <- 90 #1000, 120
y.t <- seq(-y.lim,y.lim,20)

# Set colors 
cc <- c('#1e7b59','#46d29f') # 3 up and 3 down from original

AI.s <- S(AI.a,3,6,x.t,y.t,cc) + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)

AI.s 

# Set colors 
cc <- c('#243a5b','#4975b6')

SA.s <- S(SA.a,3,6,x.t,y.t,cc) + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)


SA.s

# Set colors 
cc <- c('#4d4d00','#cccc00')

SI.s <- S(SI.a,3,6,x.t,y.t,cc) + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)


SI.s

SP <- ggarrange(TS.s,AI.s, SA.s, SI.s, ncol=2, nrow=2, align="h")
quartz()
SP

#ggsave(paste(dataset,"_scatter.pdf",sep=""), plot=SP,width = 8, height = 6, units = "in", device = 'pdf',bg = "transparent")



# Color disagreements ------------------------------------------ Color disagreements ---------------------------------------------------------------------------------------------

# Graph Function------------------------------------------------------------------------------------------------------------------------------------------------------

SL <- function(df,x.val,y.val,x.tic,y.tic,cc,x.y){
  scat <- ggplot(df, aes(x=x.val,y=y.val)) + 
    geom_point(alpha=0.75, aes(color=as.factor(df[[5]])), size=0.75) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=8, color="black"),
      text = element_text(size=10),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5),
      legend.background = element_rect(colour = "transparent",fill = "transparent")) +
    coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x="dGLS",y="ln(BF)", color="Metric Disagreement") + 
    scale_color_manual(values=cc) 
  return(scat)
}

df <- TS.a %>% 
  mutate(diff = case_when(between(TS.g,-0.5,0.5) & TS.b > 5 ~ "BF strong dGLS neutral",
                          between(TS.g,-0.5,0.5) & TS.b < -5 ~ "BF strong dGLS neutral",
                          between(TS.b,-5,5) & TS.g > 0.5 ~ "dGLS strong BF neutral",
                          between(TS.b,-5,5) & TS.g < -0.5 ~ "dGLS strong BF neutral",
                          TS.b < -5 & TS.g > 0.5 ~ "opposite strong support",
                          TS.b > 5 & TS.g < -0.5 ~ "opposite strong support",
                          TS.b > 5 & TS.g > 0.5 ~ "agree",
                          TS.b < -5 & TS.g < -0.5 ~ "agree",
                          between(TS.b,-5,5) & between(TS.g,-0.5,0.5) ~ "agree"))


max(c(abs(df$TS.g),abs(df$TS.b)))
lower <- 20 # 45
upper <- 20
x.t <- seq(-lower,upper,5) # 15
y.t <- seq(-lower,upper,5)


x.v1 <- df$TS.g
y.v1 <- df$TS.b
cc1 <- c('grey80','orange3','slateblue','black')
#quartz()
SP <- SL(df,x.v1,y.v1,x.t,y.t,cc1,x.y1) + 
  geom_vline(xintercept=c(0.5,-0.5),color=c("black"), linetype="dashed", size=0.25) +
  geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.25)

SP

#ggsave(paste(dataset,"_scatter_disagreeTS.pdf",sep=""), plot=SP,width = 6, height = 4, units = "in", device = 'pdf',bg = "transparent")
ggsave(paste(dataset,"_scatter_disagreeTSzoom.pdf",sep=""), plot=SP,width = 5, height = 3, units = "in", device = 'pdf',bg = "transparent")

