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


# Colors 
# https://www.w3schools.com/colors/colors_picker.asp
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00

# Add column to name type of support 
mLGL1 <- mutate(mLGL1, dataset = case_when(TvS != 'a' ~ "Singhal")) %>%
  select('Locus','TvS','dataset') %>% rename(TvS.g1 =  TvS)
mLGL2 <- mutate(mLGL2, dataset = case_when(TvS != 'a' ~ "SinghalOG")) %>%
  select('Locus','TvS','dataset') %>% rename(TvS.g2 =  TvS)

mLBF1 <- mutate(mLBF1, dataset = case_when(TvS != 'a' ~ "Singhal")) %>%
  select('Locus','TvS','dataset') %>% rename(TvS.b1 =  TvS)
mLBF2 <- mutate(mLBF2, dataset = case_when(TvS != 'a' ~ "SinghalOG")) %>%
  select('Locus','TvS','dataset') %>% rename(TvS.b2 =  TvS)

mLGL <- merge(mLGL1, mLGL2, by = c("Locus"))
mLBF <- merge(mLBF1, mLBF2, by = c("Locus"))


## Scatter graph


S <- function(df,xcol,ycol,x.tic,y.tic,cc,xl,yl){
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]])) + 
    geom_point(alpha=1, color=cc, size=0.5) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=14, color="black"),
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
    labs(x=xl,y=yl)
  
  return(scat)
}

limit <- 90
tic <- seq(-limit,limit,10)
cc <- color_TP
xl <- 'Singhal - data added - dGLS'
yl <- 'Singhal - original data - dGLS'
#quartz()
g <- S(mLGL,2,4,tic,tic,cc,xl,yl) +geom_abline(color=c("black"), size=0.4, linetype="dashed") 
g

limit <- 200
tic <- seq(-limit,limit,50)
cc <- color_TP
xl <- 'Singhal - data added - 2ln(BF)'
yl <- 'Singhal - original data - 2ln(BF)'
#quartz()
b <- S(mLBF,2,4,tic,tic,cc,xl,yl) +geom_abline(color=c("black"), size=0.4, linetype="dashed") 


a <- ggarrange(g,b, ncol=1, nrow=2, align="v")
a

ggsave(paste(dataset,"_scatter_",h,".pdf",sep=""), plot=a,width = 6, height = 12, units = "in", device = 'pdf',bg = "transparent")


