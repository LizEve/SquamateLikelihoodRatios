rm(list=ls())
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library(praise)
library(svglite)
library(viridis)
praise()

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
rm(list=setdiff(ls(), c("all","all.80m","all.80mb","all.60m","all.60mb")))

# Smash for scatter 

df.m <- full_join(full_join(all.60m,all.80m,by=c('Locus', "variable", "dataSet", "diff"), suffix = c("60m","80m")), all, by=c('Locus', "variable", "dataSet", "diff"))

df.mb <- full_join(full_join(all.60mb,all.80mb,by=c('Locus', "variable", "dataSet", "diff"), suffix = c("60mb","80mb")), all, by=c('Locus', "variable", "dataSet", "diff"))


# Linear Regression Functions------------------------------------------------------------------------------------------------------------------------------------------------------

lm_eqn = function(d){
  df <- d %>% group_by(variable)
  yvals <- df$GL
  xvals <- yvals %>% group_by(variable)
  m = lm(yvals ~ xvals) #y ~ x
  eq <- substitute(italic(y) == a + b %.% italic(x)*",  "~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2), 
                        b = format(unname(coef(m)[2]), digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

mytheme <- theme_bw() + theme(panel.border = element_blank()) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=6, color="black"),
    text = element_text(size=10),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(), # get rid of major grid
    plot.title = element_text(hjust = 0.5),
    legend.background = element_rect(colour = "transparent",fill = "transparent"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank()) 


# Linear Scatter Tox -----------------------------------------------------------------------------------------------------------------------------------------------------


#magma(5)
# https://www.hexcolortool.com/#b5367e
quartz()


# dataset 
df1 <- df.m %>% filter(variable != "TvS",dataSet == "Burbrink")
scat1 <- ggplot(df1) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 

df2 <- df.m %>% filter(variable != "TvS",dataSet == "Reeder")
scat2 <- ggplot(df2) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 

df3 <- df.m %>% filter(variable != "TvS",dataSet == "Singhal")
scat3 <- ggplot(df3) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 

df4 <- df.m %>% filter(variable != "TvS",dataSet == "Streicher")
scat4 <- ggplot(df4) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet) + theme(strip.background = element_blank())

s <- ggarrange(scat1,scat2,scat3,scat4, ncol=4, nrow=1,align="hv", common.legend = TRUE, legend = "right")

fig <- annotate_figure(s,
                       bottom = text_grob("dGLS", size = 12),
                       left = text_grob("ln(BF)",  size = 12, rot = 90))
fig

ggsave("Scatter_tox_bydataset.pdf", plot=fig,width = 7, height = 4, units = "in", device = 'pdf',bg = "transparent")

# TvS -----------------------------------------------------------------------------------------------------------------------------------------------------

df <- df.m %>% filter(variable == "TvS")
df$dataSet <- factor(df$dataSet, levels = c('Burbrink','Reeder','Singhal','Streicher'))
scat <- ggplot(df) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet) + theme(strip.background = element_blank())
scat

fig <- annotate_figure(scat,
                       bottom = text_grob("dGLS", size = 12),
                       left = text_grob("ln(BF)",  size = 12, rot = 90))
fig

#ggsave("Scatter_tvs_bydataset.pdf", plot=fig,width = 7, height = 2, units = "in", device = 'pdf',bg = "transparent")


#OR

# dataset 
df1 <- df.m %>% filter(variable == "TvS",dataSet == "Burbrink")
scat1 <- ggplot(df1) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 

df2 <- df.m %>% filter(variable == "TvS",dataSet == "Reeder")
scat2 <- ggplot(df2) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 

df3 <- df.m %>% filter(variable == "TvS",dataSet == "Singhal")
scat3 <- ggplot(df3) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet)+ theme(strip.background = element_blank(),strip.text.y = element_blank()) 

df4 <- df.m %>% filter(variable == "TvS",dataSet == "Streicher")
scat4 <- ggplot(df4) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.25,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.25,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.25,se=F, color="#323189",fullrange = T) +
  mytheme +
  facet_grid(variable~dataSet) + theme(strip.background = element_blank())

s <- ggarrange(scat1,scat2,scat3,scat4, ncol=4, nrow=1,align="hv", common.legend = TRUE, legend = "right")

fig <- annotate_figure(s,
                       bottom = text_grob("dGLS", size = 12),
                       left = text_grob("ln(BF)",  size = 12, rot = 90))
fig

ggsave("Scatter_tvs_bydataset.pdf", plot=fig,width = 7, height = 2, units = "in", device = 'pdf',bg = "transparent")

