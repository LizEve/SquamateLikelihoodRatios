rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
# set dataset 

dataset <- "Streicher"

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

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Histogram ---------------------------------------------- Histogram --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

H <- function(df,xval,cc,hbin,xlab){
  ggplot(data=df, aes(x=xval)) +
    geom_histogram(bins=hbin, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1) +
    #geom_histogram(breaks=brx, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1)+ 
    #geom_histogram(binwidth = max(abs(xval))*0.01, alpha=1, position="dodge",color="grey", fill=cc, size=0.1)+ 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=10),
          axis.text = element_text(size=6, color="black"),
          text = element_text(size=8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8), 
          panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid = element_blank()) +
    labs(y="Number of Loci",x=xlab)  #+ 
  #coord_cartesian(xlim = x.tic) +
  #coord_cartesian(ylim=ytic, xlim = xtic) +
  #scale_y_continuous() +
  #scale_x_continuous(breaks = xtic) +
  # geom_vline(xintercept=lines.g,color=c("black"), linetype="dashed", size=0.5) +
  #geom_vline(xintercept=c(0),color=c("black"), linetype="dashed", size=0.2) 
}


quartz()

h.bin <- 50
h <- "TvS"
#-------GL-----------------------------------------------------------------------------------------------------------------------------------------------------

# Set input data 

x.val <- TvS$GL.s
df <- TvS
cc <- 'slateblue'
# Get max min for graph to set x axis values
max.x <- round_any(max(abs(x.val)),10,f=ceiling)
x.tic <- seq(-max.x,max.x,5)
x.tic

gl <- H(df,x.val,cc,h.bin,'dGLS values') 
gl

y.tic <- seq(0,800,200)
x.tic <- seq(-10,10,5)

gl <- H(df,x.val,cc,h.bin,'dGLS values') +
  coord_cartesian(xlim = x.tic, ylim=y.tic) + 
  scale_x_continuous(breaks = x.tic) + 
  scale_y_continuous(breaks = y.tic)
gl

#--------BF-----------------------------------------------------------------------------------------------------------------------------------------------------
rm(x.val)

# Set
x.val <- TvS$BF.s
cc <- 'orange3'

h.bin <- 50

bf <- H(df,x.val,cc,h.bin,'ln(BF) values') + 
  coord_cartesian(xlim = x.tic, ylim=y.tic) + 
  scale_x_continuous(breaks = x.tic) + 
  scale_y_continuous(breaks = y.tic)
#quartz()
bf


title <- paste(dataset," - ",h,sep="") 

g <- annotate_figure(ggarrange(gl,bf, ncol=1, nrow=2, align="v"),
                     top = text_grob(title, color = "black", size = 16))
g <- g + theme_transparent()

g




#------------------------------------------------------------------------------------------------------------------------------------------------------

# Violin ---------------------------------------------- Violin --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------


V <- function(df,xval,yval,ylablab,cc,title){
  v <- ggplot(df, aes(x=xval, y=yval, fill=xval)) + 
    geom_violin(trim=TRUE) +  
    geom_point(shape = 18,size=0.5, position = position_jitterdodge(), color='black',alpha=1)+
    scale_color_manual(values=cc) + scale_fill_manual(values=cc) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=10, color="black"),
          text = element_text(size=14),
          legend.position = "none") +
    labs(y=ylablab,x="") + 
    ggtitle(title)
  return(v)
}

color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00
cc <- c(color_AI,color_SA,color_SI,color_TP)

ylab <- "standardized values"

df.g <- all[all$supportType == 'dGLS',]
xval.g <- df.g$variable
yval.g <- df.g$value

g.v <- V(df.g,xval.g,yval.g,ylab,cc,'dGLS')
g.v

df.b <- all[all$supportType == 'BF',]
xval.b <- df.b$variable
yval.b <- df.b$value

b.v <- V(df.b,xval.b,yval.b,ylab,cc,'ln(BF)') 
b.v

ytic <- c(seq(-24,24,2))

bv <- b.v + coord_cartesian(ylim=ytic) +
  scale_y_continuous(breaks = ytic)
gv <- g.v + coord_cartesian(ylim=ytic) +
  scale_y_continuous(breaks = ytic)

v <- ggarrange(gv,bv, ncol=2, nrow=1, align="v")

v

############ By metric shows better what I want 
df.a <- all[all$variable == 'AIvSA',]
xval.a <- df.a$supportType
yval.a <- df.a$value
cc <- c('slateblue','orange3')
a.v <- V(df.a,xval.a,yval.a,ylab,cc,'AIvSA') 
a.v

V <- function(df,xval,yval,ylablab,cc,title){
  v <- ggplot(df, aes(x=xval, y=yval, fill=xval)) + 
    geom_violin(trim=TRUE) +  
    geom_point(shape = 18,size=0.5, position = position_jitterdodge(), color='black',alpha=1)+
    scale_color_manual(values=cc) + scale_fill_manual(values=cc) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=10, color="black"),
          text = element_text(size=14),
          legend.position = "none") +
    labs(y=ylablab,x="") + 
    ggtitle(title)
  return(v)
}


#------------------------------------------------------------------------------------------------------------------------------------------------------

# BoxPlot ---------------------------------------------- BoxPlot --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

B <- function(df,xval,yval,ylablab,cc,title){
  b <- ggplot(df, aes(x=xval, y=yval, fill=xval)) + 
    geom_boxplot(outlier.colour="black", outlier.shape=16,
                 outlier.size=2, notch=FALSE) +
    scale_color_manual(values=cc) + scale_fill_manual(values=cc) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=10, color="black"),
          text = element_text(size=14),
          legend.position = "none") +
    labs(y=ylablab,x="") + 
    ggtitle(title)
  return(b)
}

ylab <- "standardized values"
cc <- c('slateblue','orange3')
df.a <- all[all$variable == 'AIvSA',]
xval.a <- df.a$supportType
yval.a <- df.a$value
cc <- c('slateblue','orange3')
a.v <- B(df.a,xval.a,yval.a,ylab,cc,'AIvSA') 
a.v




#------------------------------------------------------------------------------------------------------------------------------------------------------

# Histogram ---------------------------------------------- Histogram --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

H <- function(df,xval,cc,hbin,xlab){
  ggplot(data=df, aes(x=xval)) +
    geom_histogram(bins=hbin, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1) +
    #geom_histogram(breaks=brx, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1)+ 
    #geom_histogram(binwidth = max(abs(xval))*0.01, alpha=1, position="dodge",color="grey", fill=cc, size=0.1)+ 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=10),
          axis.text = element_text(size=6, color="black"),
          text = element_text(size=8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8), 
          panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid = element_blank()) +
    labs(y="Number of Loci",x=xlab)  #+ 
  #coord_cartesian(xlim = x.tic) +
  #coord_cartesian(ylim=ytic, xlim = xtic) +
  #scale_y_continuous() +
  #scale_x_continuous(breaks = xtic) +
  # geom_vline(xintercept=lines.g,color=c("black"), linetype="dashed", size=0.5) +
  #geom_vline(xintercept=c(0),color=c("black"), linetype="dashed", size=0.2) 
}



df <- all %>% filter(variable == 'AIvSI')

hbin <- 100
ggplot(df,aes(x=value)) + 
  geom_histogram(data=subset(df,supportType == 'dGLS'),fill = 'slateblue', alpha = 0.5,bins=hbin) +
  geom_histogram(data=subset(df,supportType == 'lnBF'),fill = 'orange3', alpha = 0.5,bins=hbin) 



quartz()

h.bin <- 50
h <- "TvS"
#-------GL-----------------------------------------------------------------------------------------------------------------------------------------------------

# Set input data 
h <- "TvS"
df <- all %>% filter(variable == 'TvS' & supportType=='dGLS')
x.val <- df$value
cc <- 'slateblue'
# Get max min for graph to set x axis values
max(abs(x.val))
x.tic <- seq(-10,10,5)


gl <- H(df,x.val,cc,h.bin,'dGLS values') 
gl

y.tic <- seq(0,800,200)
x.tic <- seq(-10,10,5)

gl <- H(df,x.val,cc,h.bin,'dGLS values') +
  coord_cartesian(xlim = x.tic, ylim=y.tic) + 
  scale_x_continuous(breaks = x.tic) + 
  scale_y_continuous(breaks = y.tic)
gl

#--------BF-----------------------------------------------------------------------------------------------------------------------------------------------------


# Set
df.b <- all %>% filter(variable == 'TvS' & supportType=='lnBF')
x.val.b <- df.b$value

cc.b <- 'orange3'

h.bin <- 50

bf <- H(df.b,x.val.b,cc.b,h.bin,'ln(BF) values') + 
  coord_cartesian(xlim = x.tic, ylim=y.tic) + 
  scale_x_continuous(breaks = x.tic) + 
  scale_y_continuous(breaks = y.tic)
#quartz()
bf


title <- paste(dataset," - ",h,sep="") 

g <- annotate_figure(ggarrange(gl,bf, ncol=1, nrow=2, align="v"),
                     top = text_grob(title, color = "black", size = 16))
g <- g + theme_transparent()

g




# ------------------------------------------------------------------------------------------------------------------------------------------------------

# Linear Regression ------------------------------------------ Scatter ---------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------------------




lm_eqn = function(xvals,yvals){
  m = lm(yvals ~ xvals) #y ~ x
  eq <- substitute(italic(y) == a + b %.% italic(x)*",  "~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2), 
                        b = format(unname(coef(m)[2]), digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}



df <- all.exp %>% filter(variable == 'TvS')


max(c(abs(df$value.b),abs(df$value.g)))
lower <- 10 # 
upper <- 10
x.tic <- seq(-lower,upper,5) 
y.tic <- seq(-lower,upper,5)

d.eqn <- data.frame(
  GL = -1,
  BF = 5,
  eqn = lm_eqn(df$value.g,df$value.b))

scat <- ggplot(df, aes(x=df$value.g,y=df$value.b)) + 
  geom_point(alpha=1, size=0.5, color="springgreen4") + theme_bw() + 
  theme(panel.border = element_blank()) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(), # get rid of major grid
    plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim=y.tic,xlim = x.tic) +
  scale_y_continuous(breaks = y.tic) + 
  scale_x_continuous(breaks = x.tic) +
  labs(x='dGLS',y='ln(BF)') +
  geom_smooth(method=lm ,color="black", formula = y ~ x,size=0.25,se=TRUE) 

scat + geom_text(data = d.eqn, aes(x = GL, y = BF, label = eqn),color="black", size=4, parse = TRUE) 

#----VIOLIN color outliers--------------------------------------------------------------------------------------------------------------------------------------------------

V <- function(df,xval,yval,ylablab,cc,title){
  v <- ggplot(df, aes(x=xval, y=yval, fill=xval)) + 
    geom_violin(trim=TRUE) + scale_fill_manual(values=cc) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=10, color="black"),
          text = element_text(size=14)) +
    labs(y=ylablab,x="") + 
    ggtitle(title)
  return(v)
}

ylab <- "standardized values"
h <- "AIvSA"
df.a <- a %>% filter(variable == h)
xval.a <- df.a$supportType
yval.a <- df.a$value
cc <- c('slateblue','orange3') 
a.v <- V(df.a,xval.a,yval.a,ylab,cc,h) 
a.v + geom_point(shape = 19,size=1, position = position_jitter(), aes(color=df.a$z.outlier),alpha=0.75)+
  scale_color_manual(values=cc1) 



