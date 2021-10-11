rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)


# Upload ml and bf calcs 

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh2021April.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/")

# remove individual datasets 
rm(list=setdiff(ls(), c("all")))

# Remove columns about rank 
#all <- all %>% select(-c(drank.g,drank.b,drank.diff.n))


# First I am going to collect and squish data

mdatB <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/Burbrink_iqdataNA.txt",sep=""),
                   header=TRUE) %>% drop_na() %>% mutate(dataSet = case_when(Locus != 'a' ~ "Burbrink"))


mdatR <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/Reeder_iqdataNA.txt",sep=""),
                   header=TRUE) %>% drop_na() %>% mutate(dataSet = case_when(Locus != 'a' ~ "Reeder"))

mdatSi <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/Singhal_iqdataNA.txt",sep=""),
                   header=TRUE) %>% drop_na() %>% mutate(dataSet = case_when(Locus != 'a' ~ "Singhal"))


mdatSt <- read.table(paste("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/Streicher_iqdataNA.txt",sep=""),
                   header=TRUE) %>% drop_na() %>% mutate(dataSet = case_when(Locus != 'a' ~ "Streicher"))

mdat <- rbind(rbind(rbind(mdatB,mdatR),mdatSi),mdatSt)


# Merge datasets 

# only using the ones that are in bf gl dataset 
df.all <- left_join(all,mdat,by=c("Locus","dataSet"))

# Add diff in norm vals column 

df.all$diff.norm.g.b <- (as.numeric(df.all$z.GL) - df.all$z.BF)
df.all$diff.norm.a <- abs(df.all$diff.norm.g.b)

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

r2place = function(xvals,yvals){
  xx <- range(xvals)
  yy <- range(yvals)
  x.eq <- max(xx)-mean(xx)/2
  y.eq <- min(yy)+(max(yy)-min(yy))/2
  out <- list("xax"=x.eq,"yax"=y.eq)
}
# https://www.rdocumentation.org/packages/ggpubr/versions/0.4.0/topics/ggscatter
# https://rdrr.io/cran/ggpubr/man/stat_cor.html
# https://stackoverflow.com/questions/60143052/how-to-add-r2-for-each-facet-of-ggplot-in-r

# Data for graphs------------------------------------------------------------------------------------------------------------------------------------------------------
my.pallet <- c( "#1F7ED0","#090C68", "#5C8714", "#A95307")

df1 <- df.all %>% filter(dataSet == "Burbrink") 
df2 <- df.all %>% filter(dataSet == "Reeder") 
df3 <- df.all %>% filter(dataSet == "Singhal") 
df4 <- df.all %>% filter(dataSet == "Streicher") 

# 1 Graph nBP (columns) v norm diff------------------------------------------------------------------------------------------------------------------------------------------------------

st1 <- ggscatter(df1, y = "diff.norm.a", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 4, label.x = 400, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

st2 <- ggscatter(df2, y = "diff.norm.a", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 2, label.x = 400, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st3 <- ggscatter(df3, y = "diff.norm.a", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 5, label.x = 100, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st4 <- ggscatter(df4, y = "diff.norm.a", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label = ..rr.label..),
           label.y = 9, label.x = 100, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

#Locus len and num taxa by gl, bf, norm difference, norm diff in rank
#quartz()


st <- ggarrange(st1,st2,st3,st4, ncol=4, nrow=1,align="hv", legend = "none")

st

fig <- annotate_figure(st,
                       bottom = text_grob("Base pairs per locus",  size = 10),
                       right = text_grob("AIvSA                              AIvSI                              SAvSI                              TvS",  size = 8, rot = 270),
                       left = text_grob("Absolute difference of normalized values \n dGLS - ln(BF)",  size = 10, rot = 90))
fig
ggsave("normdiff_bp.pdf", plot=fig,width = 6.938, height = 5.5, units = "in", device = 'pdf',bg = "transparent")


# 2 BF v norm diff------------------------------------------------------------------------------------------------------------------------------------------------------

st1 <- ggscatter(df1, y = "BF", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 30, label.x = 400, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 
st1
st2 <- ggscatter(df2, y = "BF", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 20, label.x = 400, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st3 <- ggscatter(df3, y = "BF", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 20, label.x = 100, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st4 <- ggscatter(df4, y = "BF", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label = ..rr.label..),
           label.y = 50, label.x = 100, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

#Locus len and num taxa by gl, bf, norm difference, norm diff in rank
#quartz()


st <- ggarrange(st1,st2,st3,st4, ncol=4, nrow=1,align="hv", legend = "none")

st

fig <- annotate_figure(st,
                       bottom = text_grob("Base pairs per locus",  size = 10),
                       right = text_grob("AIvSA                              AIvSI                              SAvSI                              TvS",  size = 8, rot = 270),
                       left = text_grob("ln(BF)",  size = 10, rot = 90))
fig
#ggsave("normdiff_BF.pdf", plot=fig,width = 6.938, height = 5.5, units = "in", device = 'pdf',bg = "transparent")




# 2 dGLS v norm diff------------------------------------------------------------------------------------------------------------------------------------------------------

st1 <- ggscatter(df1, y = "GL", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 10, label.x = 400, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

st2 <- ggscatter(df2, y = "GL", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 6, label.x = 400, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st3 <- ggscatter(df3, y = "GL", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label =  ..rr.label..),
           label.y = 15, label.x = 100, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st4 <- ggscatter(df4, y = "GL", x = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet,scales = "free") + 
  stat_cor(aes(color = variable,label = ..rr.label..),
           label.y = 20, label.x = 100, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

#Locus len and num taxa by gl, bf, norm difference, norm diff in rank
#quartz()


st <- ggarrange(st1,st2,st3,st4, ncol=4, nrow=1,align="hv", legend = "none")

st

fig <- annotate_figure(st,
                       bottom = text_grob("Base pairs per locus",  size = 10),
                       right = text_grob("AIvSA                              AIvSI                              SAvSI                              TvS",  size = 8, rot = 270),
                       left = text_grob("dGLS",  size = 10, rot = 90))
fig
ggsave("normdiff_dGLS.pdf", plot=fig,width = 6.938, height = 5.5, units = "in", device = 'pdf',bg = "transparent")














-----------------------------------------------------------------------------------------------------------------------------------------------------
###### OLD ##### -----------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------


# 1 Graph sequences(ntaxa) v norm diff------------------------------------------------------------------------------------------------------------------------------------------------------

st1 <- ggscatter(df1, x = "diff.norm.a", y = "Sequences",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3, label.y = 215, size=2, r.accuracy = 0.001, p.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

st1
st2 <- ggscatter(df2, x = "diff.norm.a", y = "Sequences",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 1.5, label.y = 95, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st3 <- ggscatter(df3, x = "diff.norm.a", y = "Sequences",
               color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 3.75, label.y = 20, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st4 <- ggscatter(df4, x = "diff.norm.a", y = "Sequences",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 5.8, label.y = 23, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

#Locus len and num taxa by gl, bf, norm difference, norm diff in rank
quartz()


st <- ggarrange(st1,st2,st3,st4, ncol=4, nrow=1,align="hv", legend = "none")

#st

fig <- annotate_figure(st,
                       left = text_grob("Number of tips",  size = 10, rot = 90),
                       right = text_grob("AIvSA                              AIvSI                              SAvSI                              TvS",  size = 8, rot = 270),
                       bottom = text_grob("Absolute difference of normalized values \n dGLS - ln(BF)",  size = 10))
fig
ggsave("normdiff_taxa.pdf", plot=fig,width = 6.938, height = 5.5, units = "in", device = 'pdf',bg = "transparent")


# 2 Graph Cols(bp) v norm diff------------------------------------------------------------------------------------------------------------------------------------------------------

sl1 <- ggscatter(df1, x = "diff.norm.a", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 3, label.y = 300, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl2 <- ggscatter(df2, x = "diff.norm.a", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 1.5, label.y = 300, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl3 <- ggscatter(df3, x = "diff.norm.a", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 3.75, label.y = 100, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl4 <- ggscatter(df4, x = "diff.norm.a", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 5.75, label.y = 100, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

#Locus len and num taxa by gl, bf, norm difference, norm diff in rank
#quartz()


sl <- ggarrange(sl1,sl2,sl3,sl4, ncol=4, nrow=1,align="hv", legend = "none")
#sl

fig2 <- annotate_figure(sl,
                       left = text_grob("Basepairs per locus",  size = 10, rot = 90),
                       right = text_grob("AIvSA                              AIvSI                              SAvSI                              TvS",  size = 8, rot = 270),
                       bottom = text_grob("Absolute difference of normalized values \n dGLS - ln(BF)",  size = 10))
fig2
ggsave("normdiff_bp.pdf", plot=fig2,width = 6.938, height = 5.5, units = "in", device = 'pdf',bg = "transparent")


# 3 Graph sequences(ntaxa) v rank diff------------------------------------------------------------------------------------------------------------------------------------------------------

st1 <- ggscatter(df1, x = "drank.diff.n", y = "Sequences",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 50, label.y = 200, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st2 <- ggscatter(df2, x = "drank.diff.n", y = "Sequences",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 40, label.y = 90, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st3 <- ggscatter(df3, x = "drank.diff.n", y = "Sequences",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 55, label.y = 15, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


st4 <- ggscatter(df4, x = "drank.diff.n", y = "Sequences",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 60, label.y = 19, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 

#quartz()

st <- ggarrange(st1,st2,st3,st4, ncol=4, nrow=1,align="hv", legend = "none")

#st

fig <- annotate_figure(st,
                       left = text_grob("Number of tips",  size = 10, rot = 90),
                       right = text_grob("AIvSA                              AIvSI                              SAvSI                              TvS",  size = 8, rot = 270),
                       bottom = text_grob("Normalized absolute difference in rank \n dGLS - ln(BF)",  size = 10))
fig
ggsave("rankdiff_taxa.pdf", plot=fig,width = 6.938, height = 5.5, units = "in", device = 'pdf',bg = "transparent")


# 4 Graph Cols(bp) v norm diff------------------------------------------------------------------------------------------------------------------------------------------------------

sl1 <- ggscatter(df1, x = "BF", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 50, label.y = 215, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 
sl1

sl2 <- ggscatter(df2, x = "drank.diff.n", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 40, label.y = 95, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl3 <- ggscatter(df3, x = "drank.diff.n", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 55, label.y = 20, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl4 <- ggscatter(df4, x = "drank.diff.n", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 60, label.y = 22, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl <- ggarrange(sl1,sl2,sl3,sl4, ncol=4, nrow=1,align="hv", legend = "none")

#sl
fig2 <- annotate_figure(sl,
                        left = text_grob("Basepairs per locus",  size = 10, rot = 90),
                        right = text_grob("AIvSA                              AIvSI                              SAvSI                              TvS",  size = 8, rot = 270),
                        bottom = text_grob("Normalized bsolute difference in rank \n dGLS - ln(BF)",  size = 10))
#fig2
#ggsave("rankdiff_bp.pdf", plot=fig2,width = 6.938, height = 5.5, units = "in", device = 'pdf',bg = "transparent")


# 5 Graph Cols(bp) v total BF------------------------------------------------------------------------------------------------------------------------------------------------------

sl1 <- ggscatter(df1, x = "drank.diff.n", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 50, label.y = 215, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl2 <- ggscatter(df2, x = "drank.diff.n", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 40, label.y = 95, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl3 <- ggscatter(df3, x = "drank.diff.n", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 55, label.y = 20, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl4 <- ggscatter(df4, x = "drank.diff.n", y = "Columns",
                 color = "variable", alpha=0.5,size=0.5,palette = my.pallet, 
                 add = "reg.line", fullrange = TRUE) + mytheme +
  facet_grid(variable~dataSet) + 
  stat_cor(aes(color = variable,label = ..rr.label..),label.x = 60, label.y = 22, size=2, r.accuracy = 0.001)+
  theme(strip.background = element_blank(),strip.text.y = element_blank()) 


sl <- ggarrange(sl1,sl2,sl3,sl4, ncol=4, nrow=1,align="hv", legend = "none")

#sl
fig2 <- annotate_figure(sl,
                        left = text_grob("Basepairs per locus",  size = 10, rot = 90),
                        right = text_grob("AIvSA                              AIvSI                              SAvSI                              TvS",  size = 8, rot = 270),
                        bottom = text_grob("Normalized bsolute difference in rank \n dGLS - ln(BF)",  size = 10))
#fig2
#ggsave("rankdiff_bp.pdf", plot=fig2,width = 6.938, height = 5.5, units = "in", device = 'pdf',bg = "transparent")



