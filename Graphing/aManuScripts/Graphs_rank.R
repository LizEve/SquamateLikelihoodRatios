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
load("Calcs_Smoosh.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/")


# Rank ------------------------------------------------------ Rank Calcs  --------------------------------------------------------------------------------------------------------
# rank vs dense rank - dense leaves no gaps when ties occur 
# dense rank makes the point the best - not sure how/if to normalize difference in rank 
# Currently dividing rank by total number of ranks to get relative position and comparing that.

# Nov 23rd - Tested single dataset/hypothesis and it matched with grouping by datsets and hypothesis. 

comp <- all %>% dplyr::group_by(dataSet,variable) %>% 
  mutate(per.g = ntile(desc(GL),100)) %>% 
  mutate(drank.g = dense_rank(desc(GL))) %>%
  mutate(rank.g = rank(desc(GL))) %>%
  mutate(per.b = ntile(desc(BF),100)) %>% 
  mutate(drank.b = dense_rank(desc(BF))) %>%
  mutate(rank.b = rank(desc(BF))) %>%
  mutate(per.diff = abs(per.g-per.b)) %>%
  mutate(drank.diff = abs(drank.g - drank.b)) %>%
  mutate(drank.diff.n = abs((drank.g/max(drank.g))-(drank.b/max(drank.b)))*100)%>%
  mutate(rank.diff = abs(rank.g - rank.b)) %>%
  mutate(rank.diff.n = abs((rank.g/max(rank.g))-(rank.b/max(rank.b)))*100)

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

df1 <- comp %>% filter(variable!="TvS", dataSet == "Burbrink")
r1 <- SC(df1,"Difference in rank","drank",c(0,100),".diff.n", "rank",1) +
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 
  #facet_wrap(~df1$variable) + mytheme 

df2 <- comp %>% filter(variable!="TvS", dataSet == "Reeder")
r2 <- SC(df2,"Difference in rank","drank",c(0,100),".diff.n", "rank",1) + 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 
  #facet_wrap(~df2$variable) + mytheme 

df3 <- comp %>% filter(variable!="TvS", dataSet == "Singhal")
r3 <- SC(df3,"Difference in rank","drank",c(0,100),".diff.n", "rank",0.5) + 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 
  #facet_wrap(~df3$variable) + mytheme 

df4 <- comp %>% filter(variable!="TvS", dataSet == "Streicher")
r4 <- SC(df4,"Difference in rank","drank",c(0,100),".diff.n", "rank",0.5) + 
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank())
  #facet_wrap(~df4$variable) + mytheme 


s <- ggarrange(r1,r2,r3,r4, ncol=4, nrow=1,align="hv", common.legend = TRUE, legend = "top")

fig <- annotate_figure(s,
                       bottom = text_grob("dGLS rank", size = 12),
                       left = text_grob("ln(BF) rank",  size = 12, rot = 90))
fig
ggsave("rank_scatter_tox_bydataset.pdf", plot=fig,width = 7, height = 4, units = "in", device = 'pdf',bg = "transparent")

# TvS -----------------------------------------------------------------------------------------------------------------------------------------------------


df1 <- comp %>% filter(variable=="TvS",dataSet == "Burbrink")
r1 <- SC(df1,"Difference in rank","drank",c(0,100),".diff.n", "rank",1)  + 
  theme(strip.background = element_blank(),strip.text.y = element_blank()) + 
  facet_grid(variable~dataSet)

df2 <- comp %>% filter(variable=="TvS",dataSet == "Reeder")
r2 <- SC(df2,"Difference in rank","drank",c(0,100),".diff.n", "rank",1) + 
  theme(strip.background = element_blank(),strip.text.y = element_blank())+ 
  facet_grid(variable~dataSet)

df3 <- comp %>% filter(variable=="TvS",dataSet == "Singhal")
r3 <- SC(df3,"Difference in rank","drank",c(0,100),".diff.n", "rank",0.5)  + 
  theme(strip.background = element_blank(),strip.text.y = element_blank()) + 
  facet_grid(variable~dataSet)

df4 <- comp %>% filter(variable=="TvS",dataSet == "Streicher")
r4 <- SC(df4,"Difference in rank","drank",c(0,100),".diff.n", "rank",0.5)  + 
  theme(strip.background = element_blank())+ 
  facet_grid(variable~dataSet)


s <- ggarrange(r1,r2,r3,r4, ncol=4, nrow=1,align="hv", legend = "none")

fig <- annotate_figure(s,
                       bottom = text_grob("dGLS rank", size = 12),
                       left = text_grob("ln(BF) rank",  size = 12, rot = 90))
fig
ggsave("rank_scatter_TvS_bydataset.pdf", plot=fig,width = 7, height = 2, units = "in", device = 'pdf',bg = "transparent")


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


