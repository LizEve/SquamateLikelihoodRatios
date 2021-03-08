

rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(readxl)



setwd("/Users/ChatNoir/Projects/Squam/Graphs/")

# Read in data
mid <- read_excel("/Users/chatnoir/Projects/Squam/Graphs/Tables/Corr_rand_pvals_mmb.xlsx")

totalLoci <- max(mid$Loci)
# Edit dataframes 
midD <- mid %>% filter(Sample == "80m" | Sample == "60m") %>% mutate(Sample = sub('sig', '', Sample))
midD$Subsample <- paste(midD$Hypothesis,midD$Sample,sep = " - ")
# factor 
midD$Subsample <- factor(midD$Subsample, levels = c("AIvSA - 80m", "AIvSI - 80m", "SAvSI - 80m", "AIvSA - 60m",  "AIvSI - 60m", "SAvSI - 60m"))
# add proportion of loci in subsample 
midD$propLoci <- midD$Loci/totalLoci


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
          text = element_text(size=10),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6), 
          panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid = element_blank()) +
    labs(y="Proportion of Loci",x=xlab) + 
    coord_cartesian(xlim = c(0,1))
}


#quartz()

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Streicher ---------------------------------------------- Streicher --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsampleStreicher")
results <- results %>% filter(Sample != "all")

# Remove non subsampled values 
mid <- midD %>% filter(Sample != "all") %>% filter(dataset == "Streicher")

# Set up numbers - TOX
df <- results %>% filter(Loci == 2000) %>% filter(Hypothesis != "TvS") 

m <- mid %>% filter(Hypothesis != "TvS") 

x.val <- df$corr.eff 
# 95% pval 
h05 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))


lines <- m$corr.eff 

# Graph 
h.bin <- 30
cc <- 'grey'
max.x <- round_any(max(abs(x.val)),10,f=ceiling)

gl <- H(df,x.val,cc,h.bin,'Correlation coefficient') + 
  geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title <- "2000 loci"

g <- gl + geom_point(m, mapping = aes(x=corr.eff, y=0, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title) 

# Set up numbers - TVS

dft2 <- results %>% filter(Loci == 300) %>% filter(Hypothesis == "TvS") 
mt2 <- mid %>% filter(Hypothesis == "TvS") 
mt2$Subsample <- paste(mt2$Hypothesis,mt2$Sample,sep = " - ")
# factor 
mt2$Subsample <- factor(mt2$Subsample, levels = c("TvS - 80m", "TvS - 60m"))
x.val2 <- dft2$corr.eff 
# 95% pval 
h05 <- unlist(dft2 %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(dft2 %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))

lines2 <- mt2$corr.eff 

# Graph 
h.bin <-30
cc <- 'grey'
max.x <- round_any(max(abs(x.val2)),10,f=ceiling)

gl2 <- H(dft2,x.val2,cc,h.bin,'Correlation coefficient') + 
  geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2) + theme(axis.title.y=element_blank())

title2 <- "300 loci"

g2 <- gl2 + geom_point(mt2, mapping = aes(x=corr.eff, y=0, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title2) 
g2

gg <- ggarrange(g,g2, ncol=2, nrow=1,align="hv", common.legend = F, legend = "right")

gg

# redoign graphs march 2021, mb values were misentered as m values in the original spreadsheet i used to make this graph 
#ggsave("CorrcoeffNull_Streicher.pdf", plot=gg,width = 7, height = 2, units = "in", device = 'pdf',bg = "transparent")


#------------------------------------------------------------------------------------------------------------------------------------------------------

# Singhal ---------------------------------------------- Singhal --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsampleSinghal")
results <- results %>% filter(Sample != "all")

# Remove non subsampled values 
mid <- midD %>% filter(Sample != "all") %>% filter(dataset == "Singhal")

# Set up numbers - TOX
df <- results %>% filter(Loci >= 3000) %>% filter(Hypothesis != "TvS") 
#d <- null2k4k %>% filter(Hypothesis != "TvS") 
m <- mid %>% filter(Hypothesis != "TvS") 


x.val <- df$corr.eff 
# 95% pval 
h05 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))


lines <- m$corr.eff 

# Graph 
h.bin <- 10
cc <- 'grey'
max.x <- round_any(max(abs(x.val)),10,f=ceiling)

gl <- H(df,x.val,cc,h.bin,'Correlation coefficient') + 
  geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title <- "3000 loci"

g <- gl + geom_point(m, mapping = aes(x=corr.eff, y=0, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title) 

# Set up numbers - TVS

dft2 <- results %>% filter(Loci == 1000) %>% filter(Hypothesis == "TvS") 
mt2 <- mid %>% filter(Hypothesis == "TvS") 
mt2$Subsample <- paste(mt2$Hypothesis,mt2$Sample,sep = " - ")
# factor 
mt2$Subsample <- factor(mt2$Subsample, levels = c("TvS - 80m", "TvS - 60m"))
x.val2 <- dft2$corr.eff 
# 95% pval 
h05 <- unlist(dft2 %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(dft2 %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))

lines2 <- mt2$corr.eff 


# Graph 
h.bin <-40
cc <- 'grey'
max.x <- round_any(max(abs(x.val2)),10,f=ceiling)

gl2 <- H(dft2,x.val2,cc,h.bin,'Correlation coefficient') + 
  geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title2 <- "1000 loci"

g2 <- gl2 + geom_point(mt2, mapping = aes(x=corr.eff, y=0, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title2) + theme(axis.title.y=element_blank())

gg <- ggarrange(g,g2, ncol=2, nrow=1,align="hv", common.legend = F, legend = "right")
gg
#ggsave("CorrcoeffNull_Singhal.pdf", plot=gg,width = 7, height = 2, units = "in", device = 'pdf',bg = "transparent")

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Reeder ---------------------------------------------- Reeder --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsampleReeder")
results <- results %>% filter(Sample != "all")

# Remove non subsampled values 
mid <- midD %>% filter(Sample != "all") %>% filter(dataset == "Reeder")

# Set up numbers - TOX
df <- results %>% filter(Loci >= 30) %>% filter(Hypothesis != "TvS") 
#d <- null2k4k %>% filter(Hypothesis != "TvS") 
m <- mid %>% filter(Hypothesis != "TvS") 

x.val <- df$corr.eff 
# 95% pval 
h05 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))


lines <- m$corr.eff 

# Graph 
h.bin <- 10
cc <- 'grey'
max.x <- round_any(max(abs(x.val)),10,f=ceiling)

gl <- H(df,x.val,cc,h.bin,'Correlation coefficient') + geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title <- "30 loci"

g <- gl + geom_point(m, mapping = aes(x=corr.eff, y=0, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title) 

# Set up numbers - TVS

dft2 <- results %>% filter(Loci == 10) %>% filter(Hypothesis == "TvS") 
mt2 <- mid %>% filter(Hypothesis == "TvS") 
mt2$Subsample <- paste(mt2$Hypothesis,mt2$Sample,sep = " - ")
# factor 
mt2$Subsample <- factor(mt2$Subsample, levels = c("TvS - 80m", "TvS - 60m"))
x.val2 <- dft2$corr.eff 
# 95% pval 
h05 <- unlist(dft2 %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(dft2 %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))

lines2 <- mt2$corr.eff 

# Graph 
h.bin <-40
cc <- 'grey'
max.x <- round_any(max(abs(x.val2)),10,f=ceiling)

gl2 <- H(dft2,x.val2,cc,h.bin,'Correlation coefficient') + 
  geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title2 <- "10 loci"

g2 <- gl2 + geom_point(mt2, mapping = aes(x=corr.eff, y=0, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title2) + theme(axis.title.y=element_blank())
g2

gg <- ggarrange(g,g2, ncol=2, nrow=1,align="hv", common.legend = F, legend = "right")
gg
# these graphs are the same as before but it doesnt make sense with the p vals
#ggsave("CorrcoeffNull_Reeder.pdf", plot=gg,width = 7, height = 2, units = "in", device = 'pdf',bg = "transparent")


#------------------------------------------------------------------------------------------------------------------------------------------------------

# Burbrink ---------------------------------------------- Burbrink --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/chatnoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/corsubsampleBurbrink")
results <- results %>% filter(Sample != "all")

# Remove non subsampled values 
mid <- midD %>% filter(Sample != "all") %>% filter(dataset == "Burbrink")

# Set up numbers - TOX
df <- results %>% filter(Loci == 200) %>% filter(Hypothesis == "AIvSA") 
#d <- null2k4k %>% filter(Hypothesis != "TvS") 
m <- mid %>% filter(Hypothesis == "AIvSA") 

x.val <- df$corr.eff 
# 95% pval 
h05 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(df %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))

lines <- m$corr.eff 

# Graph 
h.bin <- 10
cc <- 'grey'
max.x <- round_any(max(abs(x.val)),10,f=ceiling)

gl <- H(df,x.val,cc,h.bin,'Correlation coefficient') + geom_vline(xintercept=h05,color=c("black"), linetype="dashed", size=0.2)

title <- "200 loci"

g <- gl + geom_point(m, mapping = aes(x=corr.eff, y=0, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title) 
g

# Set up numbers - TVS

dft2 <- results %>% filter(Loci >= 100) %>% filter(Hypothesis == "TvS") 
mt2 <- mid %>% filter(Hypothesis == "TvS") 
mt2$Subsample <- paste(mt2$Hypothesis,mt2$Sample,sep = " - ")
# factor 
mt2$Subsample <- factor(mt2$Subsample, levels = c("TvS - 80m", "TvS - 60m"))
x.val2 <- dft2$corr.eff 
# 95% pval 
h05 <- unlist(dft2 %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.975,0.0275))))
h01 <- unlist(dft2 %>% summarise(Bin_q.95 = quantile(corr.eff, c(0.995,0.005))))

lines2 <- mt2$corr.eff 


# Graph 
h.bin <-40
cc <- 'grey'
max.x <- round_any(max(abs(x.val2)),10,f=ceiling)

gl2 <- H(dft2,x.val2,cc,h.bin,'Correlation coefficient') + 
  geom_vline(xintercept=h01,color=c("black"), linetype="dashed", size=0.2)

title2 <- "100 loci"

g2 <- gl2 + geom_point(mt2, mapping = aes(x=corr.eff, y=0, color=Subsample), position = position_dodge(width = 0.01)) + 
  ggtitle(title2) + theme(axis.title.y=element_blank())

gg <- ggarrange(g,g2, ncol=2, nrow=1,align="hv", common.legend = F, legend = "right")
gg

#ggsave("CorrcoeffNull_Burbrink.pdf", plot=gg,width = 7, height = 2, units = "in", device = 'pdf',bg = "transparent")
