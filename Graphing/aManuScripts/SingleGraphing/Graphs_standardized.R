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
# set dataset 

dataset <- "Singhal"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))


# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS")) %>%
  mutate_if(is.numeric,round,2)
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "lnBF")) %>%
  mutate_if(is.numeric,round,2)

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10
mLBF[7:10] <- mLBF[7:10]/2


# Z score
z.BF <- mLBF %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

B <- melt(z.BF,id=c("Locus","supportType"))

z.GL <- mLGL %>% select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>%
  mutate_at(c("AIvSA","AIvSI","SAvSI","TvS"), ~ scale(.,center=TRUE, scale=TRUE))

G <- melt(z.GL,id=c("Locus","supportType"))

# Smoosh together again
all <- bind_rows(B,G) 

all.exp <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b"))

# Violin ---------------------------------------------- Violin --------------------------------------------------------------------------------------------------------

V_old <- function(df,xval,yval,ylablab,cc,title,sp,s){
  v <- ggplot(df, aes(x=xval, y=yval, fill=xval)) + 
    geom_violin(trim=TRUE) +  
    geom_point(shape = sp,size=s, position = position_jitterdodge(), color="grey2",alpha=0.5)+
    #geom_point(shape = 18,size=0.75, position = position_jitterdodge(), aes(fill=xval, color=xval),alpha=0.5)+
    scale_color_manual(values=cc) + scale_fill_manual(values=cc) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=12, color="black"),
          text = element_text(size=14),
          legend.position = "none") +
    labs(y=ylablab,x="") + 
    ggtitle(title)
  return(v)
}

S <- function(df,xval,yval,ylablab,cc,title,sp,s){
  v <- ggplot(df, aes(x=xval, y=yval, fill=xval, color=xval)) + 
    geom_point(shape = sp,size=s, position = position_jitter(seed=1), color="grey1",alpha=0.7)+
    #geom_point(shape = sp,size=s, position = position_jitterdodge(seed=1), color="grey1",alpha=0.7)+
    #geom_violin(trim=TRUE,alpha=0.6,position=position_dodge(width=0.5)) +  
    scale_color_manual(values=cc) + scale_fill_manual(values=cc) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=12, color="black"),
          text = element_text(size=14),
          legend.position = "none") +
    labs(y=ylablab,x="") + 
    ggtitle(title)
  return(v)
}
V <- function(df,xval,yval,ylablab,cc,title,sp,s){
  v <- ggplot(df, aes(x=xval, y=yval, fill=xval, color=xval)) + 
    #geom_point(shape = sp,size=s, position = position_jitter(seed=1), color="grey1",alpha=0.7)+
    geom_point(shape = sp,size=s, position = position_jitterdodge(seed=1,jitter.width = .5), 
               color="grey3",alpha=0.3)+
    geom_violin(trim=TRUE,alpha=0.6,position = position_dodge()) +  
    scale_color_manual(values=cc) + scale_fill_manual(values=cc) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=16, color="black"),
          text = element_text(size=16),
          legend.position = "none",panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA)) +
    labs(y=ylablab,x="") + 
    ggtitle(title)
  return(v)
}

#V_old(df.a,xval.a,yval.a,ylab,cc,h,19,0.1) 
############ By metric shows better what I want 
ylab <- "standardized values"
h <- "AIvSI"
df.a <- all %>% filter(variable == h)
#df.a <- all %>% filter(variable==h, value != 0)
xval.a <- df.a$supportType
yval.a <- df.a$value
cc <- c('slateblue','orange3') 
#a.v <- V(df.a,xval.a,yval.a,ylab,cc,h,19,2) 
a.v <- V(df.a,xval.a,yval.a,ylab,cc,h,19,0.3) 
a.v
quartz()
yt <- seq(-3,3,1) 
a.z<- V(df.a,xval.a,yval.a,ylab,cc,h,19,0.3) + coord_cartesian(ylim=c(-3,3))+
  scale_y_continuous(breaks = yt) + ggtitle("")
a.z

v <- ggarrange(a.v,a.z, ncol=2, nrow=1,align="v")

v

#ggsave(paste(dataset,h,"violin.pdf",sep="_"), plot=v,width = 10, height = 6, units = "in", device = 'pdf',bg = "transparent")
# single graph 
#ggsave(paste(dataset,h,"violin.pdf",sep="_"), plot=a.v,width = 5, height = 6, units = "in", device = 'pdf',bg = "transparent")
ggsave(paste(dataset,h,"violinSinglePPT.pdf",sep="_"), plot=a.z,width = 3, height = 4, units = "in", device = 'pdf',bg = "transparent")


# Rank ------------------------------------------------------ Rank Calcs  --------------------------------------------------------------------------------------------------------
# rank vs dense rank - dense leaves no gaps when ties occur 
# dense rank makes the point the best - not sure how/if to normalize difference in rank 
# Currently dividing rank by total number of ranks to get relative position and comparing that.

# Oct 3rd, grouping isnt working so using only subest
h <- "AIvSA"
#df <- comp %>% filter(variable==h, value.g != 0,  value.b != 0)
all.exp.df <- all.exp %>% filter(variable==h)

comp <- all.exp.df %>% dplyr::group_by(variable) %>% 
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

SC <- function(df,x.tic,y.tic,legtit,l,c.lim,dif,lab,s){
  # defining columns here to avoid repeating this over and over
  xcol <- match(paste(l,".g",sep = ''),names(df))
  ycol <- match(paste(l,".b",sep = ''),names(df))
  CC <- match(paste(l,dif,sep = ''),names(df))
  
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]])) + 
    geom_point(alpha=1, aes(color=df[[CC]]), size=s) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=16, color="black"),
      text = element_text(size=16),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5),
      legend.background = element_rect(colour = "transparent",fill = "transparent"),
      legend.title = element_text(size=16),
      legend.position="top") +
    #coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x=paste('dGLS ',lab, sep=''),y=paste('BF ',lab,sep=''),color=legtit) + 
    guides(colour = guide_legend(override.aes = list(size=2))) +
    scale_color_viridis(limits=c.lim,direction = -1,option='D')
  #scale_color_gradientn(colors=heat.colors(2),limits=c(0,100))
  #scale_color_gradient2(low="blue", high="red",midpoint=mid)
  
  return(scat)
}


quartz()
h <- "AIvSA"
#df <- comp %>% filter(variable==h, value.g != 0,  value.b != 0)
df <- comp %>% filter(variable==h)
max(df$drank.g) # x
max(df$drank.b) # y

r <- SC(df,seq(0,800,200),seq(0,1400,200),
   "Difference in rank","drank",c(0,100),
   ".diff.n", "rank",2) + theme(plot.margin = unit(c(0.1,1,0.1,0.1), "cm"))
r 

p <- SC(df,seq(0,100,20),seq(0,100,20),
   "Difference in \n percentile","per",c(0,100),
   ".diff", "percentile",2)
p

#SC(df,seq(0,4000,500),seq(0,4000,500),"difference \n in rank","rank",c(0,100),".diff.n") # this makes a butterfly

s <- ggarrange(r,p, ncol=1, nrow=2,align="h")

s

ggsave(paste(dataset,h,"rank_scatter.pdf",sep="_"), plot=s,width = 8, height = 12, units = "in", device = 'pdf',bg = "transparent")
# single graph
ggsave(paste(dataset,h,"rank_scatterPPT.pdf",sep="_"), plot=r,width = 6, height = 4, units = "in", device = 'pdf',bg = "transparent")


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


