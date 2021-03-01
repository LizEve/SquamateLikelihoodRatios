rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# Read in data and rename data
d <- function(dataset){
  load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))
  mLBF.x <- mLBF %>% mutate(dataset = case_when(TvS != 'a' ~ dataset)) %>% select("Locus","TvS","dataset")
  mLGL.x <- mLGL %>% mutate(dataset = case_when(TvS != 'a' ~ dataset)) %>% select("Locus","TvS","dataset")
  mL <- merge(mLGL.x,mLBF.x,by = c("Locus","dataset")) %>% rename(TvS.g=TvS.x) %>% rename(TvS.b=TvS.y)
  rm(mLBF,mLGL)
  return(list(mLGL.x,mLBF.x,mL))
}

# Read in data and rename data
dataset <- "SinghalOG"
setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))
mLBF.sio <- d(dataset)[[1]]
mLGL.sio <- d(dataset)[[2]]
mL.sio <- d(dataset)[[3]]

# Read in data and rename data
dataset <- "Singhal"
mLBF.si <- d(dataset)[[1]]
mLGL.si <- d(dataset)[[2]]
mL.si <- d(dataset)[[3]]


# Read in data and rename data
dataset <- "Streicher"
mLBF.st <- d(dataset)[[1]]
mLGL.st <- d(dataset)[[2]]
mL.st <- d(dataset)[[3]]


# Read in data and rename data
dataset <- "Reeder"
mLBF.re <- d(dataset)[[1]]
mLGL.re <- d(dataset)[[2]]
mL.re <- d(dataset)[[3]]


# Read in data and rename data
dataset <- "Burbrink"
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,"_BFonly.RData", sep=''))
mLBF.bu <- mLBF %>% mutate(dataset = case_when(TvS != 'a' ~ dataset)) %>% select("Locus","TvS","dataset")
mL.bu <- mLBF.bu %>% mutate(dataset = case_when(TvS != 'a' ~ dataset)) %>% rename(TvS.b=TvS)
rm(mLBF,mLGL)


# Rename loci to all match using singhal as template
mL.bu <- data.frame(lapply(mL.bu, function(x) {gsub("T212_", "AHE-", x)}))
mL.re <- data.frame(lapply(mL.re, function(x) {gsub("Reeder_DNA_", "gene-", x)}))




## Using one of the singhal datasets, match other datasets 
# Add all these datsets together and rename columns 
rm(x,y,z,f,h,j,n)
x <- left_join(mL.si,mL.bu,by="Locus") %>% 
  rename(TvS.b.si = TvS.b.x) %>% 
  rename(TvS.g.si = TvS.g) %>% 
  rename(TvS.b.bu = TvS.b.y) 

y <- left_join(x,mL.st,by=c("Locus"))%>% 
  rename(TvS.b.st = TvS.b) %>%
  rename(TvS.g.st = TvS.g) 

z <- left_join(y,mL.re,by=c("Locus"))%>% 
  rename(TvS.b.re = TvS.b) %>%
  rename(TvS.g.re = TvS.g) 

# Merge columns in a really stupid way 
f <- unite(z, TvS.b, c(TvS.b.bu, TvS.b.st, TvS.b.re), remove=TRUE)
h <- unite(f, TvS.g, c(TvS.g.st, TvS.g.re), remove=TRUE)
j <- unite(h, dataset.y, c(dataset.y,dataset.x.x, dataset.y.y), remove=TRUE)

n <- data.frame(lapply(j, function(x) {gsub("NA|NA_|_NA|_NA_", "", x)}))


## Scatter graph


S <- function(df,xcol,ycol,x.tic,y.tic,cc,xl,yl){
  xval <- as.numeric(as.character(df[[xcol]]))
  yval <- as.numeric(as.character(df[[ycol]]))
  print(xval)
  scat <- ggplot(df, aes(x=xval,y=yval)) + 
    geom_point(alpha=0.5, aes(color=df[[cc]]), size=2) + theme_bw() + theme(panel.border = element_blank()) +
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


limit <- 200
tic <- seq(-limit,limit,50)
xl <- 'Other Datasets - 2ln(BF)'
#yl <- 'Singhal - original data - 2ln(BF)'
yl <- 'Singhal - added data - 2ln(BF)'
#quartz()
df.n <- n[!(n$dataset.y == ""), ]
b <- S(df.n,4,5,tic,tic,6,xl,yl) +geom_abline(color=c("black"), size=0.4, linetype="dashed") 
b


limit <- 60
tic <- seq(-limit,limit,10)
xl <- 'Other Datasets - dGLS'
#yl <- 'Singhal - original data - dGLS'
yl <- 'Singhal - added data - dGLS'
#quartz()
g <- S(df.n,3,7,tic,tic,6,xl,yl) +geom_abline(color=c("black"), size=0.4, linetype="dashed") 
g


a <- ggarrange(g,b, ncol=1, nrow=2, align="v")
a

ggsave(paste("Singhal_scatter_compLoci.pdf",sep=""), plot=a,width = 9, height = 9, units = "in", device = 'pdf',bg = "transparent")


