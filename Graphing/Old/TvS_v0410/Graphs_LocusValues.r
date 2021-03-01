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

# Rename 
mL <- mLGL
ytxt <- "dGLS"
yint <- c(0.5,-0.5)


mL <- mLBF
ytxt <- "2ln(BF)"
yint <- c(10,-10)


# Names "TvS", "AIvSA", "AIvSI", "SAvAI", "SAvSI", "SIvAI", "SIvSA", "TvS_support"

color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "#C2DF23FF"


quartz()
# Set colors 
color_h0 <- color_TP
color_h1 <- color_S

graph_general <- ggplot(data=mL) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  aes(x=reorder(Locus,-TvS,sum),y=TvS, color = TvS < 0, fill = TvS < 0) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=16, color="black"),
        panel.background = element_blank(),
        text = element_text(size=20),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c(color_h0,color_h1)) + 
  scale_fill_manual(values=c(color_h0,color_h1)) 

# Customize per dataset

# only need to adjust the max limits of the graph based on the largest value
max(abs(min(mL$TvS)),abs(max(mL$TvS)))
limit <- 150
#ytic <- c(seq(-limit,limit,10),0.5,-0.5)
#ytic <- seq(-limit,limit,10)
#ytic <- c(seq(-limit,limit,50),10,-10)
ytic <- c(seq(-limit,limit,20))


graph_custom <- graph_general + 
  coord_cartesian(ylim=ytic) +
  scale_y_continuous(breaks = ytic) + 
  ggtitle(paste(dataset,"\nToxPoly vs. Sclero",sep="")) + 
  labs(y=ytxt,x='') +
  geom_hline(yintercept=yint,color=c("black"), linetype="dashed", size=0.5) 
graph_custom

# Write to file 
# change axis text to 24 element texts to 30 
ggsave(paste(dataset,"_BF_Bar.pdf",sep=""), plot=graph_custom,width = 9.5, height = 6, units = "in", device = 'pdf',bg = "transparent")

ggsave(paste(dataset,"_dGLS_Bar.pdf",sep=""), plot=graph_custom,width = 9.5, height = 6, units = "in", device = 'pdf',bg = "transparent")
