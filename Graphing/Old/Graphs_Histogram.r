rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# Colors 
# https://www.w3schools.com/colors/colors_picker.asp
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00
# set dataset 

#dataset <- "SinghalOG"
dataset <- "Streicher"


setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS")) %>%
  mutate_if(is.numeric,round,2)
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))



# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10 - not including average calcs
names(mLBF[7:10])
mLBF[7:10] <- mLBF[7:10]/2
# for singhal OG 
#names(mLBF[4])
#mLBF[4] <- mLBF[4]/2

# Get total number of loci 
loci <- length(mLGL$Locus)



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


#quartz()

h.bin <- 50
h <- "TvS"
#-------GL-----------------------------------------------------------------------------------------------------------------------------------------------------

# Set input data 

x.val <- mLGL$TvS
df <- mLGL
cc <- 'slateblue'
# Get max min for graph to set x axis values
max.x <- round_any(max(abs(x.val)),10,f=ceiling)
x.tic <- seq(-max.x,max.x,round(max.x/5,digits=-1))
max(abs(x.val))
x.tic

gl <- H(df,x.val,cc,h.bin,'dGLS values') +
  coord_cartesian(xlim = x.tic) + 
  scale_x_continuous(breaks = x.tic) 

#quartz()
gl

y.tic <- seq(0,800,200)

gl <- H(df,x.val,cc,h.bin,'dGLS values') +
  coord_cartesian(xlim = x.tic, ylim=y.tic) + 
  scale_x_continuous(breaks = x.tic) + 
  scale_y_continuous(breaks = y.tic)
gl

#--------BF-----------------------------------------------------------------------------------------------------------------------------------------------------
rm(x.val,x.tic)

# Set
x.val <- mLBF$TvS
df <- mLBF
cc <- 'orange3'
# Get max min for graph to set x axis values
max.x <- round_any(max(abs(x.val)),10,f=ceiling)
x.tic <- seq(-max.x,max.x,round(max.x/5,digits=-1))
x.tic <- seq(-max.x,max.x,10)
max(abs(x.val))
x.tic
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

#ggsave(paste(dataset,"_histo_",h,".pdf",sep=""), plot=g,width = 4, height = 3, units = "in", device = 'pdf',bg = "transparent")

#---Tx------------------------------------------------------------------------------------------------------------------------------------------------------------



Htx <- function(x.val){
  # Get max min for graph to set x axis values
  max.x <- round_any(max(abs(x.val)),10,f=ceiling)
  x.tic <- seq(-max.x,max.x,round(max.x/2,digits=0))
  return(x.tic)
  }
Htx(mLGL$AIvSA)
# Singhal has outlier that messes up the graphs - AHE-L95, AHE-L8, AHE-L354
#nope <- c('AHE-L95', 'AHE-L8', 'AHE-L354')
#mLBF <- mLBF %>% filter(!Locus %in% nope)
#mLGL <- mLGL %>% filter(!Locus %in% nope)

#----------TxGL-----------------------------------------------------------------------------------------------------------------------------------------------------
nclass.FD(mLGL$AIvSA)
quartz()
# Set input data 
hbin <- 50
y.tic <- seq(0,4000,1000)
cc <- 'slateblue'
g1 <- H(mLGL,mLGL$AIvSA,cc,hbin,'dGLS') + 
  coord_cartesian(xlim = seq(-60,30, by=10), ylim=y.tic) + 
  #scale_x_continuous(breaks=Htx(mLGL$AIvSA)) +
  scale_x_continuous(breaks=seq(-60,30, by=10)) +
  scale_y_continuous(breaks = y.tic) #+ ggtitle('AIvSA') 
g1
#x <- seq(-0.5,0.5,0.05)
#g1 + coord_cartesian(xlim = x, ylim=y.tic) + 
#  scale_x_continuous(breaks=x)
g2 <- H(mLGL,mLGL$AIvSI,cc,hbin,'dGLS') + ggtitle('AIvSI') + 
  coord_cartesian(xlim = Htx(mLGL$AIvSI), ylim=y.tic) + 
  scale_x_continuous(breaks=Htx(mLGL$AIvSI)) + 
  scale_y_continuous(breaks = y.tic)

g3 <- H(mLGL,mLGL$SAvSI,cc,hbin,'dGLS') + ggtitle('SAvSI') + 
  coord_cartesian(xlim = Htx(mLGL$SAvSI), ylim=y.tic) + 
  scale_x_continuous(breaks=Htx(mLGL$SAvSI)) + 
  scale_y_continuous(breaks = y.tic)


#-------TxBF-----------------------------------------------------------------------------------------------------------------------------------------------------
cc <- 'orange3'

b1 <- H(mLBF,mLBF$AIvSA,cc,hbin,'ln(BF)')  + 
  coord_cartesian(xlim = seq(-60,30, by=10), ylim=y.tic) + 
  scale_x_continuous(breaks=seq(-60,30, by=10)) +
  #scale_x_continuous(breaks=Htx(mLBF$AIvSA)) + 
  scale_y_continuous(breaks = y.tic) #+ ggtitle('AIvSA')
b1
b2 <- H(mLBF,mLBF$AIvSI,cc,hbin,'ln(BF)') + ggtitle('AIvSI') + 
  coord_cartesian(xlim = Htx(mLBF$AIvSI), ylim=y.tic) + 
  scale_x_continuous(breaks=Htx(mLBF$AIvSI)) + 
  scale_y_continuous(breaks = y.tic)

b3 <- H(mLBF,mLBF$SAvSI,cc,hbin,'ln(BF)') + ggtitle('SAvSI') + 
  coord_cartesian(xlim = Htx(mLBF$SAvSI), ylim=y.tic) + 
  scale_x_continuous(breaks=Htx(mLBF$SAvSI)) + 
  scale_y_continuous(breaks = y.tic)

title <- dataset 

#g <- annotate_figure(ggarrange(g1,g2,g3,b1,b2,b3, ncol=3, nrow=2, align="v"),
#                     top = text_grob(title, color = "black", size = 16))
#g

gv <- g1 + geom_vline(xintercept=c(-0.5,0.5),color=c("black"))
gv
bv <- b1 + geom_vline(xintercept=c(-5,5),color=c("black"))
bv


gv <- g1 + geom_hline(yintercept=3170,color=c("black"),linetype="dashed",size=0.2) 
bv <- b1 + geom_hline(yintercept=3170,color=c("black"),linetype="dashed",size=0.2) +
  geom_hline(yintercept=2490,color=c("black"),size=0.2)
bv

v <- ggarrange(gv,bv, ncol=1, nrow=2, align="v")
v

h <- "aivsa"
ggsave(paste(dataset,"_histo_",h,".pdf",sep=""), plot=v,width = 3, height = 3, units = "in", device = 'pdf',bg = "transparent")



### old 

## Get max number of loci around 0
#brx <- pretty(range(x.val), n = nclass.FD(x.val), min.n = 1)
#binwidth <- brx[2]-brx[1]
#max.y <- round_any(sum(x.val > 0 & x.val < binwidth),50,f=ceiling)
#y.tic <- seq(0,max.y,round(max.y/10,digits=-1))
#sum(x.val > 0 & x.val < binwidth)
#y.tic
#H <- function(df,xval,cc,xtic,brx,xlab){
#  ggplot(data=df, aes(x=xval)) +
#    geom_histogram(bins=100, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1) +
#    #geom_histogram(breaks=brx, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1)+ 
#    #geom_histogram(binwidth = max(abs(xval))*0.01, alpha=1, position="dodge",color="grey", fill=cc, size=0.1)+ 
#    theme_classic() + 
#    theme(plot.title = element_text(hjust = 0.5, size=16),
#          axis.text = element_text(size=8, color="black"),
#          text = element_text(size=12),
#          legend.title = element_text(size = 10),
#          legend.text = element_text(size = 10)) +
#    labs(y="Number of Loci",x=xlab)  + 
#    #coord_cartesian(xlim = xtic) +
#    #coord_cartesian(ylim=ytic, xlim = xtic) +
#    #scale_y_continuous() +
#    #scale_x_continuous(breaks = xtic) +
#    # geom_vline(xintercept=lines.g,color=c("black"), linetype="dashed", size=0.5) +
#    geom_vline(xintercept=c(0),color=c("black"), linetype="dashed", size=0.2) 
#}
#