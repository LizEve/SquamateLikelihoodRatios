#rm(list=ls())
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
setwd("/Users/ChatNoir/Projects/Squam/Graphs")

all.e[all.e=="BF"]<-"ln(BF)"
all[all=="BF"]<-"ln(BF)"


# Violin ---------------------------------------------- Violin --------------------------------------------------------------------------------------------------------

V <- function(df,cc,title,sp,s,a){
  v <- ggplot(df, aes(x=supportType, y=z.value, fill=supportType, color=supportType)) + 
    #geom_point(shape = sp,size=s, position = position_jitter(seed=1), color="grey1",alpha=0.7)+
    geom_point(shape = sp,size=s, position = position_jitterdodge(seed=1,jitter.width = .5), 
               color="grey3",alpha=a)+
    geom_violin(trim=TRUE,alpha=0.2,position = position_dodge()) +  
    scale_color_manual(values=cc,labels=c("dGLS","ln(BF)")) + scale_fill_manual(values=cc,labels=c("dGLS","ln(BF)")) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0, size=14),
          axis.text = element_text(size=8, color="black"),
          text = element_text(size=10),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    ggtitle(title)
  return(v)
}

# by dataset  ------------------------------------------------------ by dataset violin --------------------------------------------------------------------------------------------------------

cc <- c('slateblue','orange3') 
v <- all.e %>% filter(variable != "TvS")
v1 <- all.e %>% filter(dataSet == "Streicher" |  dataSet == "Singhal")
v2 <- all.e %>% filter(dataSet == "Burbrink" |  dataSet == "Reeder")

f1 <- V(v1,cc,"a.",19,0.3,0.1) + 
  facet_grid(dataSet ~ variable) 

f2 <- V(v2,cc,"b.",19,0.3,0.1) + 
  facet_grid(dataSet ~ variable) 

p <- ggarrange(f1,f2, ncol=1, nrow=2,align="v")

p

fig <- annotate_figure(p,
                       bottom = text_grob("Hypotheses compared", size = 14),
                       left = text_grob("Standardized values",  size = 14, rot = 90),
                       top = text_grob("Distribution of standardized values",  size = 16),
                       fig.lab = NULL)
fig

setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
#ggsave("SI_ByData_allhypoth.pdf", plot=fig,width = 7, height = 5, units = "in", device = 'pdf',bg = "transparent")


## Zoom in 
tt <- seq(-2,2,1) 
yt <- c(-2,2)

f1 <- V(v1,cc,"a.",19,0.2,0.1) + 
  facet_grid(dataSet ~ variable) + 
  coord_cartesian(ylim=yt) +
  scale_y_continuous(breaks = tt)

f2 <- V(v2,cc,"b.",19,0.2,0.1) + 
  facet_grid(dataSet ~ variable) + 
  coord_cartesian(ylim=yt) +
  scale_y_continuous(breaks = tt)

p <- ggarrange(f1,f2, ncol=1, nrow=2,align="v")

p

fig <- annotate_figure(p,
                       bottom = text_grob("Hypotheses compared", size = 14),
                       left = text_grob("Standardized values",  size = 14, rot = 90),
                       top = text_grob("Distribution of standardized values",  size = 16),
                       fig.lab = NULL)
fig

setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
#ggsave("SI_ByData_allhypoth_zoom.pdf", plot=fig,width = 7, height = 5, units = "in", device = 'pdf',bg = "transparent")

# combined data  ------------------------------------------------------ combined data violin --------------------------------------------------------------------------------------------------------

#quartz()
cc <- c('slateblue','orange3') 


v <- all.e 
f1 <- V(v,cc,"a.",19,0.3,0.1) + 
  facet_wrap(~ v$variable,ncol=4,strip.position = c("bottom")) 
# for non zoom use specific ranges for each dataset
tt <- seq(-4,4,2) 
yt <- c(-4,4)

f2 <- V(v,cc,"b.",19,0.1,0.1) + 
  facet_wrap(~ v$variable,ncol=4,strip.position = c("bottom")) + 
  coord_cartesian(ylim=yt) +
  scale_y_continuous(breaks = tt)
#f2

p <- ggarrange(f1,f2, ncol=1, nrow=2,align="v")
#p 

fig <- annotate_figure(p,
                       bottom = text_grob("Hypotheses compared", size = 14),
                       left = text_grob("Standardized values",  size = 14, rot = 90),
                       top = text_grob("Distribution of standardized values",  size = 16),
                       fig.lab = NULL)
fig

setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
#ggsave("SI_CombindedData_allhypoth.pdf", plot=fig,width = 5, height = 5, units = "in", device = 'pdf',bg = "transparent")

# single violin ------------------------------------------------------ lonely violin --------------------------------------------------------------------------------------------------------
cc <- c('slateblue','orange3') 


v <- all.e %>% filter(variable %in% c("AIvSA","TvS")) 
tt <- seq(-6,6,2) 
yt <- c(-6,6)

f1 <- V(v,cc,"",19,0.2,0.1) + 
  facet_wrap(~ Burbrink$variable,ncol=4,strip.position = c("bottom")) + 
  coord_cartesian(ylim=yt) +
  scale_y_continuous(breaks = tt)

#f1

fig <- annotate_figure(f1,
                       bottom = text_grob("Hypotheses compared", size = 14),
                       left = text_grob("Standardized values",  size = 14, rot = 90),
                       top = text_grob("Distribution of standardized values",  size = 16))
fig

setwd("/Users/ChatNoir/Projects/Squam/Graphs/")
#ggsave("CombinedData_violin_zoom.pdf", plot=fig,width = 5, height = 4, units = "in", device = 'pdf',bg = "transparent")


 # Rank ------------------------------------------------------ Rank Calcs  --------------------------------------------------------------------------------------------------------
# rank vs dense rank - dense leaves no gaps when ties occur 
# dense rank makes the point the best - not sure how/if to normalize difference in rank 
# Currently dividing rank by total number of ranks to get relative position and comparing that.

comp <- all.exp %>% group_by(variable) %>% 
  mutate(per.g = ntile(desc(value.g),100)) %>% 
  mutate(drank.g = dense_rank(desc(value.g))) %>%
  mutate(rank.g = rank(desc(value.g))) %>%
  mutate(per.b = ntile(desc(value.b),100)) %>% 
  mutate(drank.b = dense_rank(desc(value.b))) %>%
  mutate(rank.b = rank(desc(value.b))) %>%
  mutate(per.diff = abs(per.g-per.b)) %>%
  mutate(drank.diff = abs(drank.g - drank.b)) %>%
  mutate(drank.diff.n = abs((drank.g/max(drank.g))-(drank.b/max(drank.b)))*100)%>%
  mutate(rank.diff = abs(rank.g - rank.b)) %>%
  mutate(rank.diff.n = abs((rank.g/max(rank.g))-(rank.b/max(rank.b)))*100)

# Scatter Rank ---------------------------------------------- Rank Scatter --------------------------------------------------------------------------------------------------------

SC <- function(df,x.tic,y.tic,legtit,l,c.lim,dif,lab,s,tit){
  # defining columns here to avoid repeating this over and over
  xcol <- match(paste(l,".g",sep = ''),names(df))
  ycol <- match(paste(l,".b",sep = ''),names(df))
  CC <- match(paste(l,dif,sep = ''),names(df))
  
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]])) + 
    geom_point(alpha=1, aes(color=df[[CC]]), size=s) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=8, color="black"),
      text = element_text(size=10),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5, size=12),
      legend.background = element_rect(colour = "transparent",fill = "transparent"),
      legend.title = element_text(size=8)) +
    #coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x=paste('dGLS ',lab, sep=''),y=paste('ln(BF) ',lab,sep=''),color=legtit) + 
    guides(colour = guide_legend(override.aes = list(size=2))) +
    scale_color_viridis(limits=c.lim,direction = -1,option='D') + 
    ggtitle(tit)
  #scale_color_gradientn(colors=heat.colors(2),limits=c(0,100))
  #scale_color_gradient2(low="blue", high="red",midpoint=mid)
  
  return(scat)
}


#quartz()
h <- "TvS"
#df <- comp %>% filter(variable==h, value.g != 0,  value.b != 0)
df <- comp %>% filter(variable==h)
max(df$drank.g) # x
max(df$drank.b) # y
r <- SC(df,seq(0,2200,400),seq(0,2800,400),
   "Normalized \n difference \n in rank","drank",c(0,100),
   ".diff.n", "rank",1, "a. TvS")

p <- SC(df,seq(0,100,20),seq(0,100,20),
   "Difference \n in percentile","per",c(0,100),
   ".diff", "percentile",1, "c. TvS")
#p

h1 <- "AIvSA"
#df <- comp %>% filter(variable==h, value.g != 0,  value.b != 0)
df1 <- comp %>% filter(variable==h1)
df1 <- comp 
max(df1$drank.g) # x
max(df1$drank.b) # y
r1 <- SC(df1,seq(0,2200,200),seq(0,2800,400),
        "Normalized \n difference \n in rank","drank",c(0,100),
        ".diff.n", "rank",1, "b. AIvSA")
r1

p1 <- SC(df1,seq(0,100,20),seq(0,100,20),
        "Difference in \n percentile","per",c(0,100),
        ".diff", "percentile",1, "d. AIvSA")
#p1

#SC(df,seq(0,4000,500),seq(0,4000,500),"difference \n in rank","rank",c(0,100),".diff.n") # this makes a butterfly

s <- ggarrange(r,r1, ncol=2, nrow=1,align="h",common.legend = TRUE, 
               legend = "right")
s
s1 <- ggarrange(p,p1, ncol=2, nrow=1,align="h",common.legend = TRUE, 
                legend = "right")

#s1


ss <- ggarrange(s,s1, ncol=1, nrow=2,align="h")
#ss


figm <- annotate_figure(ss,
                        top = text_grob("Singhal",  size = 14))
figm


ggsave(paste("Fig4_",dataset,h,"rank_scatter.pdf",sep="_"), plot=ss,width = 8, height = 8, units = "in", device = 'pdf',bg = "transparent")
# single graph
#ggsave(paste(dataset,h,"rank_scatter_nozero.pdf",sep="_"), plot=r,width = 8, height = 6, units = "in", device = 'pdf',bg = "transparent")

#---------test
scat <- ggplot(df, aes(x=df$drank.diff.n,y=df$value.g)) + 
  geom_point(alpha=1, aes(color=df$value.b), size=1) + theme_bw() + theme(panel.border = element_blank()) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(), # get rid of major grid
    plot.title = element_text(hjust = 0.5),
    legend.background = element_rect(colour = "transparent",fill = "transparent"),
    legend.title = element_text(size=12)) +
  labs(x='n rank diff',y="gl val",color="bf val") +
  scale_color_viridis(option='D')
scat
scat <- ggplot(df, aes(x=df$drank.diff.n,y=df$drank.g)) + 
  geom_point(alpha=1, aes(color=df$drank.b), size=1) + theme_bw() + theme(panel.border = element_blank()) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(), # get rid of major grid
    plot.title = element_text(hjust = 0.5),
    legend.background = element_rect(colour = "transparent",fill = "transparent"),
    legend.title = element_text(size=12)) +
  labs(x='n rank diff',y="gl rank",color="bf rank") +
  scale_color_viridis(option='D')
scat
#------------
