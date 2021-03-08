rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
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



# Graph Function------------------------------------------------------------------------------------------------------------------------------------------------------

SL <- function(df,x.val,y.val,cc,colorby,lowerx,upperx,ticx,lowery,uppery,ticy){
  x.tic <- seq(-lowerx,upperx,ticx) 
  y.tic <- seq(-lowery,uppery,ticy)
  #x.tic <- c(-lowerx,upperx,ticx) 
  #y.tic <- c(-lowery,uppery,ticy)
  scat <- ggplot(df, aes(x=x.val,y=y.val)) + 
    geom_point(alpha=0.75, aes(color=as.factor(.data[[colorby]])), size=.75) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=12, color="black"),
      text = element_text(size=14),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.background = element_rect(colour = "transparent",fill = "transparent")) +
    #coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    #labs(x="dGLS",y="ln(BF)", color="Metric Disagreement") + 
    scale_color_manual(values=cc) 
  return(scat)
}

# Linear Regression Functions------------------------------------------------------------------------------------------------------------------------------------------------------

lm_eqn = function(xvals,yvals){
  m = lm(yvals ~ xvals) #y ~ x
  eq <- substitute(italic(y) == a + b %.% italic(x)*",  "~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2), 
                        b = format(unname(coef(m)[2]), digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

LR <- function(dd,h,cc,colorby,lowerx,upperx,ticx,lowery,uppery,ticy,xe,ye,ds){
  df <- dd %>% filter(variable == h)
  x.val <- df$GL
  y.val <- df$BF
  d.eqn <- data.frame(
    GL = xe,
    BF = ye,
    eqn = lm_eqn(df$GL,df$BF))
  
  lr <- SL(df,x.val,y.val,cc,colorby,lowerx,upperx,ticx,lowery,uppery,ticy) +
    geom_smooth(method=lm ,color="black", formula = y ~ x,size=0.25,se=TRUE) + 
    #geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.25) +
    #geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.25) + 
    theme(legend.position = "none") + ggtitle(ds) + 
    geom_text(data = d.eqn, aes(x = mean(range(x.val)), y = BF, label = eqn),color="black", 
              size=4, parse = TRUE) +
    geom_vline(xintercept=0,color=c("black"), linetype="dashed", size=0.25) + 
    geom_hline(yintercept=0,color=c("black"), linetype="dashed", size=0.25)
  return(lr) 
}



# Linear Scatter TvS -----------------------------------------------------------------------------------------------------------------------------------------------------

dataset <- "All"

tt <- paste(dataset,"Toxicofera v Scleroglossa",sep = " - ")
h <- 'TvS'

df <- all %>% filter(variable == h)
xval <- df$GL
yval <- df$BF
cc1 <- c(  "midnightblue","slateblue1", "#B63679FF", "#FB8861FF", "#FCFDBFFF")
#magma(5)
quartz()
# x y
a.lr <- LR(all,h,cc1,'dataSet',80,80,20,150,150,50,0,160,"All")
a.lr

a.lr.80m <- LR(all.80m,h,cc1,'dataSet',60,60,10,60,50,10,0,50,"80m") 
a.lr.80m

a.lr.80mb <- LR(all.80mb,h,cc1,'dataSet',60,60,1,60,50,5,0,10,"80mb") 
a.lr.80mb

a.lr.60m <- LR(all.60m,h,cc1,'dataSet',60,60,10,60,50,10,0,40,"60m") 
a.lr.60m

a.lr.60mb <- LR(all.60mb,h,cc1,'dataSet',60,60,1,60,50,5,0,5,"60mb") 
a.lr.60mb




lr <- ggarrange(a.lr.80m,a.lr.80mb, a.lr.60m,a.lr.60mb,a.lr,  ncol=2, nrow=3,align="h", common.legend = TRUE, legend='bottom')
lr


fig <- annotate_figure(lr,
                       bottom = text_grob("dGLS", size = 16),
                       left = text_grob("ln(BF)",  size = 16, rot = 90),
                       top = text_grob(tt,  size = 18))
fig


ggsave(paste("Alldata",h,"scatter_lrm.pdf",sep="_"), plot=fig,width = 8, height = 8, units = "in", device = 'pdf',bg = "transparent")



# Linear Scatter Tox -----------------------------------------------------------------------------------------------------------------------------------------------------
LR <- function(dd,h,cc,colorby,lowerx,upperx,ticx,lowery,uppery,ticy,xe,ye,ds){
  df <- dd 
  x.val <- df$GL
  y.val <- df$BF
  d.eqn <- data.frame(
    GL = xe,
    BF = ye,
    eqn = lm_eqn(df$GL,df$BF))
  
  lr <- SL(df,x.val,y.val,cc,colorby,lowerx,upperx,ticx,lowery,uppery,ticy) +
    geom_smooth(method=lm ,color="black", formula = y ~ x,size=0.25,se=TRUE) + 
    #geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.25) +
    #geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.25) + 
    theme(legend.position = "none") + ggtitle(ds) + 
    geom_text(data = d.eqn, aes(x = mean(range(x.val)), y = BF, label = eqn),color="black", 
              size=4, parse = TRUE) +
    geom_vline(xintercept=0,color=c("black"), linetype="dashed", size=0.25) + 
    geom_hline(yintercept=0,color=c("black"), linetype="dashed", size=0.25)
  return(lr) 
}

dataset <- "all"
#tt <- paste(dataset,"AI v SA",sep = " - ")
#h <- 'AIvSA'
#tt <- paste(dataset,"AI v SI",sep = " - ")
#h <- 'AIvSI'
#tt <- paste(dataset,"SA v SI",sep = " - ")
#h <- 'SAvSI'
#df <- all %>% filter(variable == h)

df <- all %>% filter(variable != "TvS")
h <- "inTox"

xval <- df$GL
yval <- df$BF

#quartz()
cc1 <- c( "slateblue1", "midnightblue", "#B63679FF", "#FB8861FF", "#FCFDBFFF")
a.lr <- LR(all,h,cc1,'dataSet',60,60,20,160,160,20,0,140,"")
#a.lr 
 
a.lr.80m <- LR(all.80m,h,cc1,'dataSet',50,50,5,60,60,10,0,40,"") 
#a.lr.80m

a.lr.50m <- LR(all.50m,h,cc1,'dataSet',40,40,2,60,60,5,0,40,"") 
#a.lr.50m

# named 2x in paper cuz its 2x threshold vut its +1 threshold above threshold.
a.lr.2x <- LR(all.1x,h,cc1,'dataSet',40,40,2,60,60,5,0,40,"2x") 
#a.lr.2x

a.lr.1x <- LR(all.0x,h,cc1,'dataSet',40,40,2,60,60,5,0,40,"1x") 
#a.lr.1x


lr <- ggarrange(a.lr,a.lr.80m,a.lr.2x,a.lr.1x,  ncol=2, nrow=2,align="h")
#lr

#tt <- "all data, all tox hyp"
fig <- annotate_figure(lr,
                       bottom = text_grob("dGLS", size = 16),
                       left = text_grob("ln(BF)",  size = 16, rot = 90),
                       top = text_grob(tt,  size = 18))
fig
ggsave(paste("Alldata",h,"scatter_lrm.pdf",sep="_"), plot=fig,width = 8, height = 8, units = "in", device = 'pdf',bg = "transparent")

