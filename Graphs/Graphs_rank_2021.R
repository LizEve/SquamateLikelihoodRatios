rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library("viridis")

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh2021April.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/")



# Scatter Rank ---------------------------------------------- Rank Scatter --------------------------------------------------------------------------------------------------------

SC <- function(df,legtit,l,c.lim,dif,lab,s){
  # defining columns here to avoid repeating this over and over
  xcol <- match(paste(l,".g",sep = ''),names(df))
  ycol <- match(paste(l,".b",sep = ''),names(df))
  CC <- match(paste(l,dif,sep = ''),names(df))
  
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]],fill=dataSet)) + 
    geom_point(alpha=1, aes(color=df[[CC]]), size=s, position = "identity") + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=6, color="black"),
      text = element_text(size=10),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      #legend.background = element_rect(colour = "transparent",fill = "transparent"),
      #legend.title = element_text(size=16),
      #legend.position="top",
      legend.position = "none",
      legend.title = element_blank()) +
    labs(color=legtit)+
    guides(colour = guide_legend(override.aes = list(size=2))) +
    scale_color_viridis(limits=c(0,100),direction = -1,option='D')
  
  return(scat)
}

# Tox scatter ------------------------------------------------------ Tox scatter -------------------------------------------------------------------------------------------------------

quartz()

df1 <- all %>% filter(variable!="TvS", dataSet == "Burbrink")
r1 <- SC(df1,"Difference in rank","drank",c(0,100),".diff.n", "rank",1) +
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 
  #facet_wrap(~df1$variable) + mytheme 

df2 <- all %>% filter(variable!="TvS", dataSet == "Reeder")
r2 <- SC(df2,"Difference in rank","drank",c(0,100),".diff.n", "rank",1) + 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 
  #facet_wrap(~df2$variable) + mytheme 

df3 <- all %>% filter(variable!="TvS", dataSet == "Singhal")
r3 <- SC(df3,"Difference in rank","drank",c(0,100),".diff.n", "rank",0.5) + 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 
  #facet_wrap(~df3$variable) + mytheme 

df4 <- all %>% filter(variable!="TvS", dataSet == "Streicher")
r4 <- SC(df4,"Difference in rank","drank",c(0,100),".diff.n", "rank",0.5) + 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank())
  #facet_wrap(~df4$variable) + mytheme 


r <- ggarrange(r1,r2,r3,r4, ncol=4, nrow=1,align="hv", common.legend = TRUE, legend = "top")

fig <- annotate_figure(r,
                       left = text_grob("ln(BF) rank",  size = 10, rot = 90))
fig
#ggsave("rank_scatter_tox_bydataset.pdf", plot=fig,width = 7, height = 4, units = "in", device = 'pdf',bg = "transparent")

# TvS -----------------------------------------------------------------------------------------------------------------------------------------------------


d1 <- all %>% filter(variable=="TvS",dataSet == "Burbrink")
s1 <- SC(d1,"Difference in rank","drank",c(0,100),".diff.n", "rank",1)  + 
  theme(strip.background = element_blank(),strip.text.y = element_blank()) + 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.x = element_blank())

d2 <- all %>% filter(variable=="TvS",dataSet == "Reeder")
s2 <- SC(d2,"Difference in rank","drank",c(0,100),".diff.n", "rank",1) + 
  theme(strip.background = element_blank(),strip.text.y = element_blank())+ 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.x = element_blank())

d3 <- all %>% filter(variable=="TvS",dataSet == "Singhal")
s3 <- SC(d3,"Difference in rank","drank",c(0,100),".diff.n", "rank",0.5)  + 
  theme(strip.background = element_blank(),strip.text.y = element_blank()) + 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.x = element_blank())

d4 <- all %>% filter(variable=="TvS",dataSet == "Streicher")
s4 <- SC(d4,"Difference in rank","drank",c(0,100),".diff.n", "rank",0.5)  + 
  theme(strip.background = element_blank())+ 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.x = element_blank())


s <- ggarrange(s1,s2,s3,s4, ncol=4, nrow=1,align="hv", legend = "none")


fig <- annotate_figure(s,
                       left = text_grob("ln(BF) rank",  size = 10, rot = 90))
fig
ggsave("rank_scatter_TvS_bydataset1.4.pdf", plot=fig,width = 7, height = 1.4, units = "in", device = 'pdf',bg = "transparent")
ggsave("rank_scatter_TvS_bydataset1.6.pdf", plot=fig,width = 7, height = 1.6, units = "in", device = 'pdf',bg = "transparent")


# Few scatter ------------------------------------------------------ lonely scatter -------------------------------------------------------------------------------------------------------

df1 <- comp %>% filter(variable=="AIvSA",dataSet == "Singhal")
r1 <- SC(df1,"Difference in rank","drank",c(0,100),".diff.n", "rank",1) +
  facet_wrap(~df1$variable,strip.position = c("bottom"))

df2 <- comp %>% filter(variable=="TvS",dataSet == "Singhal")
r2 <- SC(df2,"Difference in rank","drank",c(0,100),".diff.n", "rank",1) + 
  facet_wrap(~df2$variable,strip.position = c("bottom"))


s <- ggarrange(r1,r2, ncol=2, nrow=1,align="hv", common.legend = TRUE, legend = "right")

fig <- annotate_figure(s,
                       bottom = text_grob("dGLS rank", size = 12),
                       left = text_grob("ln(BF) rank",  size = 12, rot = 90))

ggsave("rank_scatter_Singhal.pdf", plot=fig,width = 7, height = 3, units = "in", device = 'pdf',bg = "transparent")


