rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(readxl)





setwd("/Users/ChatNoir/Projects/Squam/Graphs/")

# Read in data
mid <- read_excel("/Users/chatnoir/Projects/Squam/Graphs/Tables/AllSchmoosh.xlsx")
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsample")
totalLoci <- max(mid$Loci)
# Remove non subsampled values 
mid <- mid %>% filter(Sample != "all")
results <- results %>% filter(Sample != "all")

# Edit dataframes 
mid <- mid %>% filter(Sample == "sig80m" | Sample == "sig60m") %>% mutate(Sample = sub('sig', '', Sample))
mid$Subsample <- paste(mid$Hypothesis,mid$Sample,sep = " - ")
# factor 
mid$Subsample <- factor(mid$Subsample, levels = c("AIvSA - 80m", "AIvSI - 80m", "SAvSI - 80m", "AIvSA - 60m",  "AIvSI - 60m", "SAvSI - 60m"))
# add proportion of loci in subsample 
mid$propLoci <- mid$Loci/totalLoci


# Select null distributions for 2k-4k 
null2k4k <- results %>% filter(Loci >= 2000) 
null4k <- results %>% filter(Loci >= 4000) 
null2k <- results %>% filter(Loci >= 2000) %>% filter(Loci <= 4000) 


#------------------------------------------------------------------------------------------------------------------------------------------------------

# Pvals ---------------------------------------------- Pvals --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
h <- 'TvS'
n <- null2k %>% filter(Hypothesis == h) 
s <- '80m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
s <- '60m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
h <- 'AIvSA'
n <- null2k %>% filter(Hypothesis == h) 
s <- '80m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
s <- '60m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
h <- 'AIvSA'
n <- null2k %>% filter(Hypothesis == h) 
s <- '80m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
s <- '60m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
h <- 'AIvSI'
n <- null2k %>% filter(Hypothesis == h) 
s <- '80m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
s <- '60m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
h <- 'SAvSI'
n <- null2k %>% filter(Hypothesis == h) 
s <- '80m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s
s <- '60m'
m <- mid %>% filter(Hypothesis == h) %>% filter(Sample == s)
e <- m$corr.eff
p.s <- sum(abs(n$corr.eff) >= abs(e)) / length(n$corr.eff)
p.s

# all p = 1, ie all sim cor eff are higher than emp vals. all of them. 

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Histogram ---------------------------------------------- Histogram --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

H <- function(df,xval,cc,hbin,xlab){
  ggplot(data=df, aes(x=xval)) +
  geom_histogram(aes(y=..count../sum(..count..)),bins=hbin, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1) +
  #geom_histogram(breaks=brx, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1)+ 
  #geom_histogram(binwidth = max(abs(xval))*0.01, alpha=1, position="dodge",color="grey", fill=cc, size=0.1)+ 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=10),
        axis.text = element_text(size=6, color="black"),
        text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6), 
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid = element_blank()) +
  labs(y="Proportion of Loci",x=xlab) + 
    coord_cartesian(xlim = c(0,1))
  }


quartz()


# Set up numbers - TOX
df <- results %>% filter(Loci >= 4000) %>% filter(Hypothesis != "TvS") 
#d <- null2k4k %>% filter(Hypothesis != "TvS") 
m <- mid %>% filter(Hypothesis != "TvS") 

x.val <- df$corr.eff 
# 95% pval 
h05 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))


lines <- m$corr.eff 

# Graph 
h.bin <- 20
cc <- 'grey'
max.x <- round_any(max(abs(x.val)),10,f=ceiling)

gl <- H(df,x.val,cc,h.bin,'Correlation coefficient') + geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title <- "4000 loci"

g <- gl + geom_point(m, mapping = aes(x=corr.eff, y=p, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title) 
g


# Set up numbers - TVS
dft <- results %>% filter(Loci == 2000) %>% filter(Hypothesis == "TvS") 
mt <- mid %>% filter(Hypothesis == "TvS") 
mt$Subsample <- paste(mt$Hypothesis,mt$Sample,sep = " - ")
# factor 
mt$Subsample <- factor(mt$Subsample, levels = c("TvS - 80m", "TvS - 60m"))
x.val <- dft$corr.eff 
# 95% pval 
h05 <- unlist(dft %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(dft %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))


lines <- mt$corr.eff 

# Graph 
h.bin <-40
cc <- 'grey'
max.x <- round_any(max(abs(x.val)),10,f=ceiling)

gl2 <- H(dft,x.val,cc,h.bin,'Correlation coefficient') + geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title <- "2000 loci"

g2 <- gl2 + geom_point(mt, mapping = aes(x=corr.eff, y=p, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title) 

gg <- ggarrange(g,g2, ncol=2, nrow=1,align="hv", common.legend = F, legend = "right")


ggsave("CorrcoeffNull_all.pdf", plot=gg,width = 7, height = 2, units = "in", device = 'pdf',bg = "transparent")



#------------------------------------------------------------------------------------------------------------------------------------------------------

# two types subsample -----------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------




setwd("/Users/ChatNoir/Projects/Squam/Graphs/")

# Read in data
mid <- read_excel("/Users/chatnoir/Projects/Squam/Graphs/Tables/AllSchmoosh_mb.xlsx")
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsample")
totalLoci <- max(mid$Loci)
# Remove non subsampled values 
mid <- mid %>% filter(Sample != "all")
results <- results %>% filter(Sample != "all")

# Edit dataframes 
mid <- mid %>%  mutate(Sample = sub('quant', '', Sample))
mid$Subsample <- paste(mid$Hypothesis,mid$Sample,sep = " - ")
# factor 
mid$Subsample <- factor(mid$Subsample, levels = c("AIvSA - 80m",  "AIvSI - 80m", "SAvSI - 80m", 
                                                  "AIvSA - 80mb",  "AIvSI - 80mb", "SAvSI - 80mb", 
                                                  "AIvSA - 60m",  "AIvSI - 60m", "SAvSI - 60m",
                                                  "AIvSA - 60mb",  "AIvSI - 60mb", "SAvSI - 60mb"))
# add proportion of loci in subsample 
mid$propLoci <- mid$Loci/totalLoci


#------------------------------------------------------------------------------------------------------------------------------------------------------

# Histogram ---------------------------------------------- Histogram --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

quartz()


# Set up numbers - TOX
df <- results %>% filter(Loci >= 4000) %>% filter(Hypothesis != "TvS") 
#d <- null2k4k %>% filter(Hypothesis != "TvS") 
m <- mid %>% filter(Hypothesis != "TvS") 

x.val <- df$corr.eff 
# 95% pval 
h05 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))


lines <- m$corr.eff 

# Graph 
h.bin <- 20
cc <- 'grey'
max.x <- round_any(max(abs(x.val)),10,f=ceiling)

gl <- H(df,x.val,cc,h.bin,'Correlation coefficient') + geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title <- "Null correlation coefficient distribution for random subsample of 4k loci"

g <- gl + geom_point(m, mapping = aes(x=corr.eff, y=p, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title) 
g
#ggsave("CorrcoeffNull_tox_all4k_mb.pdf", plot=g,width = 6, height = 4, units = "in", device = 'pdf',bg = "transparent")


# Set up numbers - TVS

dft <- results %>% filter(Loci == 1000) %>% filter(Hypothesis == "TvS") 
mt <- mid %>% filter(Hypothesis == "TvS") 
mt$Subsample <- paste(mt$Hypothesis,mt$Sample,sep = " - ")
# factor 
mt$Subsample <- factor(mt$Subsample, levels = c("TvS - 80m","TvS - 80mb", "TvS - 60m", "TvS - 60mb"))
x.val <- dft$corr.eff 
# 95% pval 
h05 <- unlist(dft %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(dft %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))


lines <- mt$corr.eff 

# Graph 
h.bin <-40
cc <- 'grey'
max.x <- round_any(max(abs(x.val)),10,f=ceiling)

gl <- H(dft,x.val,cc,h.bin,'Correlation coefficient') + geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title <- "Null correlation coefficient distribution for random subsample of 1k loci"

g <- gl + geom_point(mt, mapping = aes(x=corr.eff, y=p, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title) 
g


#ggsave("CorrcoeffNull_tvs_all1k_mb.pdf", plot=g,width = 6, height = 4, units = "in", device = 'pdf',bg = "transparent")



