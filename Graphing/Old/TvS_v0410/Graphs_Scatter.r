rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# set dataset 

dataset <- "Reeder"


setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Merge datasets, then select and rename columns of interest 

x <- merge(mLBF, mLGL, by = c("Locus"))
mL <- x %>% select("Locus")

mL$BF <- x$TvS.x
mL$dGLS <- x$TvS.y
hypoth <- "TvS"
title <- "ToxPoly vs. Sclero"

#mL$BF <- x$AIvSA.x
#mL$dGLS <- x$AIvSA.y
#hypoth <- "AIvSA"
#title <- "ToxAI vs. ToxSA"
#

# Get max min for graph 

max(abs(min(mL$BF)),abs(max(mL$BF)),abs(min(mL$dGLS)),abs(max(mL$dGLS)))
limit <- 150
tic <- seq(-limit,limit,20)

# Names "TvS", "AIvSA", "AIvSI", "SAvAI", "SAvSI", "SIvAI", "SIvSA", "TvS_support"

color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "#C2DF23FF"


quartz()
# Set colors 
color_h0 <- color_AI
color_h1 <- color_S

graph_general <- ggplot(mL, aes(x=BF,y=dGLS)) + 
  geom_point(alpha=0.5, color=color_h0, size=1) + theme_bw() + theme(panel.border = element_blank()) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=16, color="black"),
    text = element_text(size=20),
    legend.position = "none", 
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(), # get rid of major grid
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim=tic,xlim = tic) +
  scale_y_continuous(breaks = tic) + 
  scale_x_continuous(breaks = tic) +
  labs(x='2ln(BF)',y='dGLS')


graph_custom <- graph_general + 
  geom_vline(xintercept=c(10,-10),color=c("gray"), linetype="dotted", size=0.5) +
  geom_hline(yintercept=c(0.5,-0.5),color=c("gray"), linetype="dotted", size=0.5) +
  #geom_abline(color=c("gray"), size=0.2, linetype="solid") + 
  ggtitle(paste(dataset,"\n",title,sep=""))

graph_custom

ggsave(paste(dataset,"_scatter_",hypoth,".pdf",sep=""), plot=graph_custom,width = 9.5, height = 6, units = "in", device = 'pdf',bg = "transparent")
