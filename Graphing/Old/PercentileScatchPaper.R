rm(list=ls())
library(plyr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)



# https://stackoverflow.com/questions/10317972/percentile-plot-with-ggplot2-bars-with-y-and-yend
# https://stackoverflow.com/questions/35096287/ggplot2-boxplot-with-geometric-mean-and-90th-and-10th-percentiles
# https://stackoverflow.com/questions/4765482/changing-whisker-definition-in-geom-boxplot/4765608#4765608
# https://stackoverflow.com/questions/3494593/shading-a-kernel-density-plot-between-two-points/4371473#4371473
# https://stackoverflow.com/questions/14863744/adding-percentile-lines-to-a-density-plot/14864404#14864404
# https://stackoverflow.com/questions/34029811/fill-different-colors-for-each-quantile-in-geom-density-of-ggplot
# https://stackoverflow.com/questions/30488389/using-dplyr-window-functions-to-calculate-percentiles


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



mLGL <- mutate(mLGL, decile_rank = ntile(mLGL$TvS,10))
mLGL <- mutate(mLGL, centile_rank = ntile(mLGL$TvS,100))
z <- mLGL %>% group_by(decile_rank) %>% dplyr::summarise(d_counts=n())



histogram colored by percentile 

limit <- round(max(abs(min(mLBF$TvS)),abs(max(mLBF$TvS))),-1)
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
quantiles <- quantile(mLBF$TvS, prob=probs) # same numbers as plyr call
mLBF$quantTvS <- factor(findInterval(mLBF$TvS,quantiles))
ggplot(mLBF, aes(TvS, fill=mLBF$quantTvS, color=mLBF$quantTvS)) + 
  geom_histogram(binwidth = limit*0.01, alpha=1, position="dodge")
quartz()

limit <- round(max(abs(min(mLGL$TvS)),abs(max(mLGL$TvS))),-1)
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
quantiles <- quantile(mLGL$TvS, prob=probs) # same numbers as plyr call
mLGL$quantTvS <- factor(findInterval(mLGL$TvS,quantiles))
ggplot(mLGL, aes(TvS, fill=mLGL$quantTvS, color=mLGL$quantTvS)) + 
  geom_histogram(binwidth = limit*0.01, alpha=1, position="dodge")









# Set title and input data 
title.g <- paste(dataset," - Tox vs Sclero",sep="") 
cc.g <- color_TP
df.g <- mLBF
xval.g <- mLBF$TvS
xlablab.g <- 'dGLS values'
per <- seq(0,1,0.1)



x <- df.g %>% summarise(TvS = list(enframe(quantile(TvS, probs=per)))) %>% unnest


dt <- data.table(x=c(1:200),y=rnorm(200))
dens <- density(dt$y)
df <- data.frame(x=dens$x, y=dens$y)
probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
quantiles <- quantile(dt$y, prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))
ggplot(df, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + scale_x_continuous(breaks=quantiles) + scale_fill_brewer(guide="none")



datSum <- ddply(df.g,.(supportType),summarise,q05 = quantile(xval.g,0.05),
                q25 = quantile(xval.g,0.25),
                q45 = quantile(xval.g,0.45),
                q55 = quantile(xval.g,0.55),
                q75 = quantile(xval.g,0.75),
                q95 = quantile(xval.g,0.95))
datSum$NameAlt <- 1:1
quartz()
ggplot(datSum,aes(ymin = NameAlt - 0.2,ymax = NameAlt + 0.2)) + 
  geom_rect(aes(xmin = q05,xmax = q95),fill = "white",colour = "black") + 
  geom_rect(aes(xmin = q25,xmax = q75),fill = 'lightblue',colour = "black") +
  geom_rect(aes(xmin = q45,xmax = q55),fill = 'blue',colour = "black") +
  scale_y_continuous(breaks = 1:1,labels = datSum$Name)

