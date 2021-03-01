rm(list=ls())
library(plyr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
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

# Calculate percentile numbers
per <- seq(0,1,0.1)
x <- mLGL %>% summarise(TvS = list(enframe(quantile(TvS, probs=per)))) %>% unnest

# Set title and input data 
title.g <- paste(dataset," - AI v SA",sep="") 
cc.g <- color_AI
df.g <- mLGL
xval.g <- mLGL$AIvSA
xlablab.g <- 'dGLS'

datSum.g <- ddply(df.g,.(supportType),summarise,
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
                q9 = quantile(xval.g,0.9))
datSum.g$NameAlt <- 1:1
datSum.g
tics.g <- seq(-0.7,0.8,0.01) # Streicher 

GL <- ggplot(datSum.g,aes(ymin = NameAlt - 0.2,ymax = NameAlt + 0.2)) + 
  geom_rect(aes(xmin = q1,xmax = q9),fill = "white",colour = "black") + 
  geom_rect(aes(xmin = q2,xmax = q8),fill = "lightblue",colour = "black") + 
  geom_rect(aes(xmin = q3,xmax = q7),fill = "dodgerblue",colour = "black") + 
  geom_rect(aes(xmin = q4,xmax = q6),fill = 'blue',colour = "black") +
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


# Set title and input data 
title <- paste(dataset," - Tox vs Sclero",sep="") 
cc <- color_TP
df <- mLBF
xval <- mLBF$TvS
xlablab <- '2ln(BF)'

datSum <- ddply(df,.(supportType),summarise,
                q1 = quantile(xval,0.1),
                q2 = quantile(xval,0.2),
                q3 = quantile(xval,0.3),
                q4 = quantile(xval,0.4),
                q5a = quantile(xval,0.49),
                q5 = quantile(xval,0.5),
                q5b = quantile(xval,0.51),
                q6 = quantile(xval,0.6),
                q7 = quantile(xval,0.7),
                q8 = quantile(xval,0.8),
                q9 = quantile(xval,0.9))
datSum$NameAlt <- 1:1
datSum
tics <- seq(-10,27,5) # Streicher 
BF <- ggplot(datSum,aes(ymin = NameAlt - 0.2,ymax = NameAlt + 0.2)) + 
  geom_rect(aes(xmin = q1,xmax = q9),fill = "white",colour = "black") + 
  geom_rect(aes(xmin = q2,xmax = q8),fill = "lightblue",colour = "black") + 
  geom_rect(aes(xmin = q3,xmax = q7),fill = "dodgerblue",colour = "black") + 
  geom_rect(aes(xmin = q4,xmax = q6),fill = 'blue',colour = "black") +
  #geom_rect(aes(xmin = q5,xmax = q5),fill = 'midnightblue',colour = "black") +
  theme_classic() + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size=16),
        axis.text = element_text(size=10, color="black"),
        text = element_text(size=14)) + 
  labs(x=xlablab,y="") +
  coord_cartesian(xlim=tics) +
  scale_x_continuous(breaks = tics)

graph_combined <- ggarrange(BF, GL, ncol=1, nrow=2, align="v")
x <- annotate_figure(graph_combined, top = text_grob(paste("Percentiles - ",dataset,sep = ""),size = 16))

ggsave(paste(dataset,"_decile.pdf",sep=""), plot=x, width = 8, height = 5, units = "in", device = 'pdf',bg = "transparent")

