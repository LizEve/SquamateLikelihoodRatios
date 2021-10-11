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
setwd("/Users/ChatNoir/Projects/Squam/Graphs")

all.bf <- all %>% select(c(Locus,variable,z.BF,dataSet)) %>% 
  mutate(supportType = case_when(Locus != 'a' ~ "ln(BF)")) %>% 
  dplyr::rename(z.value = z.BF) 

all.gl <- all %>% select(c(Locus,variable,z.GL,dataSet)) %>% 
  mutate(supportType = case_when(Locus != 'a' ~ "dGLS"))%>% 
  dplyr::rename(z.value = z.GL) 

all <- bind_rows(all.bf,all.gl)

# Violin ---------------------------------------------- Violin --------------------------------------------------------------------------------------------------------

V <- function(df,cc,sp,s,a){
  v <- ggplot(df, aes(x=supportType, y=z.value, fill=supportType, color=supportType)) + 
    #geom_point(shape = sp,size=s, position = position_jitter(seed=1), color="grey1",alpha=0.7)+
    geom_point(shape = sp,size=s, position = position_jitterdodge(seed=1,jitter.width = .5), 
               color="grey3",alpha=a)+
    geom_violin(trim=TRUE,alpha=0.2,position = position_dodge()) +  
    scale_color_manual(values=cc,labels=c("dGLS","ln(BF)")) + scale_fill_manual(values=cc,labels=c("dGLS","ln(BF)")) +
    theme_classic() + 
    theme(axis.text = element_text(size=6, color="black"),
          text = element_text(size=10),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) 
  return(v)
}

# by dataset  ------------------------------------------------------ by dataset violin --------------------------------------------------------------------------------------------------------

cc <- c('slateblue','orange3') 
v <- all %>% filter(variable != "TvS")
v1 <- all %>% filter(dataSet == "Streicher" |  dataSet == "Singhal")
v2 <- all %>% filter(dataSet == "Burbrink" |  dataSet == "Reeder")

f1 <- V(v1,cc,19,0.3,0.1) + 
  facet_grid(dataSet ~ variable) + theme(strip.background = element_blank()) + 
  theme(axis.line.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

f2 <- V(v2,cc,19,0.3,0.1) + 
  facet_grid(dataSet ~ variable) + theme(strip.background = element_blank(),
                                         strip.text.x = element_blank()) 


p <- ggarrange(f1,f2, ncol=1, nrow=2,align="v")

p

fig <- annotate_figure(p,
                       bottom = text_grob("Hypotheses compared", size = 10),
                       left = text_grob("Standardized values",  size = 10, rot = 90),
                       fig.lab = NULL)
fig

setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
ggsave("Violins.pdf", plot=fig,width = 7, height = 5, units = "in", device = 'pdf',bg = "transparent")


## Zoom in 
tt <- seq(-2,2,1) 
yt <- c(-2,2)

f1 <- V(v1,cc,19,0.2,0.1) + 
  facet_grid(dataSet ~ variable) + 
  coord_cartesian(ylim=yt) +
  scale_y_continuous(breaks = tt)+ theme(strip.background = element_blank()) + 
  theme(axis.line.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

f2 <- V(v2,cc,19,0.2,0.1) + 
  facet_grid(dataSet ~ variable) + 
  coord_cartesian(ylim=yt) +
  scale_y_continuous(breaks = tt) + theme(strip.background = element_blank(),
                                          strip.text.x = element_blank()) 

p <- ggarrange(f1,f2, ncol=1, nrow=2,align="v")

p

fig <- annotate_figure(p,
                       bottom = text_grob("Hypotheses compared", size = 10),
                       left = text_grob("Standardized values",  size = 10, rot = 90),,
                       fig.lab = NULL)
fig

setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
ggsave("Violins_zoom.pdf", plot=fig,width = 7, height = 5, units = "in", device = 'pdf',bg = "transparent")

# # combined data  ------------------------------------------------------ combined data violin --------------------------------------------------------------------------------------------------------
# 
# #quartz()
# cc <- c('slateblue','orange3') 
# 
# 
# v <- all.e 
# f1 <- V(v,cc,"a.",19,0.3,0.1) + 
#   facet_wrap(~ v$variable,ncol=4,strip.position = c("bottom")) 
# # for non zoom use specific ranges for each dataset
# tt <- seq(-4,4,2) 
# yt <- c(-4,4)
# 
# f2 <- V(v,cc,"b.",19,0.1,0.1) + 
#   facet_wrap(~ v$variable,ncol=4,strip.position = c("bottom")) + 
#   coord_cartesian(ylim=yt) +
#   scale_y_continuous(breaks = tt)
# #f2
# 
# p <- ggarrange(f1,f2, ncol=1, nrow=2,align="v")
# #p 
# 
# fig <- annotate_figure(p,
#                        bottom = text_grob("Hypotheses compared", size = 14),
#                        left = text_grob("Standardized values",  size = 14, rot = 90),
#                        top = text_grob("Distribution of standardized values",  size = 16),
#                        fig.lab = NULL)
# fig
# 
# setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
# #ggsave("SI_CombindedData_allhypoth.pdf", plot=fig,width = 5, height = 5, units = "in", device = 'pdf',bg = "transparent")
# 
# # single violin ------------------------------------------------------ lonely violin --------------------------------------------------------------------------------------------------------
# cc <- c('slateblue','orange3') 
# 
# 
# v <- all.e %>% filter(variable %in% c("AIvSA","TvS")) 
# tt <- seq(-6,6,2) 
# yt <- c(-6,6)
# 
# f1 <- V(v,cc,"",19,0.2,0.1) + 
#   facet_wrap(~ Burbrink$variable,ncol=4,strip.position = c("bottom")) + 
#   coord_cartesian(ylim=yt) +
#   scale_y_continuous(breaks = tt)
# 
# #f1
# 
# fig <- annotate_figure(f1,
#                        bottom = text_grob("Hypotheses compared", size = 14),
#                        left = text_grob("Standardized values",  size = 14, rot = 90),
#                        top = text_grob("Distribution of standardized values",  size = 16))
# fig
# 
# setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
# #ggsave("CombinedData_violin_zoom.pdf", plot=fig,width = 5, height = 4, units = "in", device = 'pdf',bg = "transparent")

