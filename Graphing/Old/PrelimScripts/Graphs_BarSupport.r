rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# set dataset 

dataset <- "Reeder"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Colors 
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "#C2DF23FF"

# Bar graph of inferred support bins 

# We need to name our vector of factors and lables so everything colors properly 
myList <- c("springgreen4","darkolivegreen3","gray55","goldenrod1","orange")
names(myList) <- c("Very Strong for Tox","Strong for Tox","Ambiguous","Strong for Sclero","Very Strong for Sclero")
myColors <- c("Very Strong for Tox","Strong for Tox","Ambiguous","Strong for Sclero","Very Strong for Sclero")
names(myColors) <- c("springgreen4","darkolivegreen3","gray55","goldenrod1","orange")


title <- "Loci binned by 2ln(BF) values"
legend_title <- "Interpretation of support"
xlabs <- c("10 and over ", "6 to 10", "6 to -6", "-6 to -10", "-10 and under")
xlablab <- '2ln(BF) values'
ytic <- c(seq(0,40,5))
df <- TvS_BFsupport
xval <- TvS_BFsupport$TvS_support
yval <- TvS_BFsupport$counts

BF_bar <- ggplot(data=df, aes(x=xval,y=yval)) + 
  geom_bar(aes(color = xval, fill = xval), 
           stat = "identity", position = position_dodge()) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text = element_text(size=10, color="black"),
        text = element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values=myList, guide=FALSE) + 
  labs(y="Number of Loci",x=xlablab)  +
  scale_y_continuous(breaks = ytic) + 
  scale_fill_manual(values=myList, 
                    name=legend_title,
                    limits = myColors) +
  scale_x_discrete(labels=xlabs,drop = FALSE) + 
  ggtitle(paste(dataset,"\n",title,sep=""))
BF_bar

title.g <- paste("Loci binned by dGLS values", sep='')
legend_title.g <- "Interpretation of support"
xlabs.g <- c("0.5 and over ", "0.1 to 0.5", "0.1 to -0.1", "-0.1 to -0.5", "-0.5 and under")
xlablab.g <- 'dGLS values'
ytic.g <- c(seq(0,45,5))
df.g <- TvS_GLsupport
xval.g <- TvS_GLsupport$TvS_support
yval.g <- TvS_GLsupport$counts

GL_bar <- ggplot(data=df.g, aes(x=xval.g,y=yval.g)) + 
  geom_bar(aes(color = xval.g, fill = xval.g), 
           stat = "identity", position = position_dodge()) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text = element_text(size=10, color="black"),
        text = element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values=myList, guide=FALSE) + 
  labs(y="Number of Loci",x=xlablab.g)  +
  scale_y_continuous(breaks = ytic.g) + 
  scale_fill_manual(values=myList, 
                    name=legend_title.g, 
                    limits = myColors) +
  scale_x_discrete(labels=xlabs.g,drop = FALSE) + 
  ggtitle(title.g)

GL_bar 

graph_combined <- ggarrange(BF, GL, ncol=1, nrow=2, legend="bottom", align="v",common.legend = TRUE)
graph_combined
# Write to file 
# change axis text to 24 element texts to 30 
ggsave(paste(dataset,"_support.pdf",sep=""), plot=graph_combined,width = 6, height = 9.5, units = "in", device = 'pdf',bg = "transparent")



