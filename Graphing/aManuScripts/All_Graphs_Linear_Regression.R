rm(list=ls())
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library(praise)
library(svglite)
library(viridis)
praise()

setwd("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles")
load("Calcs_Smoosh2021.RData")
setwd("/Users/ChatNoir/Projects/Squam/Graphs/")


# Need to redo subsets without splitting by dataset ------------------------------------------------------------------------------------------------------------------------------------------------------


# all bf within and all gl within. agnostic to other value. gl in OR bf in 
all.60m <- all %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) | 
           between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80m <- all %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) | 
           between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

# both within - gl AND bf within  
all.60mb <- all %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.2),quantile(GL, 0.8)) , 
         between(BF,quantile(BF,0.2),quantile(BF, 0.8)))

all.80mb <- all %>% group_by(variable) %>% 
  filter(between(GL,quantile(GL,0.1),quantile(GL, 0.9)) , 
         between(BF,quantile(BF,0.1),quantile(BF, 0.9)))

rm(list=setdiff(ls(), c("all","all.80m","all.80mb","all.60m","all.60mb")))

# Smash for scatter 

df.m <- full_join(full_join(all.60m,all.80m,by=c('Locus', "variable", "dataSet", "diff"), suffix = c("60m","80m")), all, by=c('Locus', "variable", "dataSet", "diff"))

df.mb <- full_join(full_join(all.60mb,all.80mb,by=c('Locus', "variable", "dataSet", "diff"), suffix = c("60mb","80mb")), all, by=c('Locus', "variable", "dataSet", "diff"))


# Linear Regression Functions------------------------------------------------------------------------------------------------------------------------------------------------------

lm_eqn = function(subset,xvals,yvals){
  m = lm(yvals ~ xvals) #y ~ x
  eq <- substitute(s*": "~~italic(r)^2~"="~r2, 
                   list(s = format(subset), 
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}



mytheme <- theme_bw() + theme(panel.border = element_blank()) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=8, color="black"),
    text = element_text(size=10),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(), # get rid of major grid
    plot.title = element_text(hjust = 0.5),
    legend.background = element_rect(colour = "transparent",fill = "transparent")) 


# Linear Scatter Tox -----------------------------------------------------------------------------------------------------------------------------------------------------


#magma(5)
# https://www.hexcolortool.com/#b5367e
quartz()

h <- 'TvS'

# dataset 
df <- df.m %>% filter(variable == h)
#h <- 'Tox'
#df <- df.m %>% filter(variable != "TvS")
x.tic <- seq(-80,80,10) 
y.tic <- seq(-150,150,20)

# position of equation on chart
#xe <- mean(range(df$GL))
xe <- -40

d.eqn <- data.frame(GL = xe, BF = 60, eqn = lm_eqn("all",df$GL,df$BF))
d.eqn.80 <- data.frame(GL = xe, BF = 50, eqn = lm_eqn("80m",df$GL80m,df$BF80m))
d.eqn.60 <- data.frame(GL = xe, BF = 40, eqn = lm_eqn("60m",df$GL60m,df$BF60m))

scat <- ggplot(df) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.5,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.5,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.5,se=F, color="#323189",fullrange = T) +
  mytheme +
  scale_y_continuous(breaks = y.tic) + 
  scale_x_continuous(breaks = x.tic)+
  labs(x="dGLS",y="ln(BF)") + 
  geom_text(data = d.eqn, aes(x = GL, y = BF, label = eqn),color="#830346", size=2, parse = TRUE, hjust = 0) + 
  geom_text(data = d.eqn.80, aes(x = GL, y = BF, label = eqn),color="#6956E6", size=2, parse = TRUE, hjust = 0) + 
  geom_text(data = d.eqn.60, aes(x = GL, y = BF, label = eqn),color="#323189", size=2, parse = TRUE, hjust = 0) +
  labs(x="dGLS",y="ln(BF)") 

scat


ggsave(paste("Combineddata",h,"scatter_lrm2.pdf",sep="_"), plot=scat,width = 4, height = 3, units = "in", device = 'pdf',bg = "transparent")



# TOX mb -----------------------------------------------------------------------------------------------------------------------------------------------------

# dataset 
df <- df.mb %>% filter(variable == h)
df <- df.mb %>% filter(variable != "TvS")
x.tic <- seq(-80,80,10) 
y.tic <- seq(-150,150,20)

# position of equation on chart
#xe <- mean(range(df$GL))
xe <- -30

d.eqn <- data.frame(GL = xe, BF = 60, eqn = lm_eqn(df$GL,df$BF))
d.eqn.80 <- data.frame(GL = xe, BF = 50, eqn = lm_eqn(df$GL80mb,df$BF80mb))
d.eqn.60 <- data.frame(GL = xe, BF = 40, eqn = lm_eqn(df$GL60mb,df$BF60mb))

scat <- ggplot(df) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80mb, BF80mb),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60mb, BF60mb),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.5,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80mb,BF80mb), method=lm,formula = y ~ x,size=0.5,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60mb,BF60mb), method=lm,formula = y ~ x,size=0.5,se=F, color="#323189",fullrange = T) +
  mytheme +
  scale_y_continuous(breaks = y.tic) + 
  scale_x_continuous(breaks = x.tic)+
  labs(x="dGLS",y="ln(BF)") + 
  geom_text(data = d.eqn, aes(x = GL, y = BF, label = eqn),color="#830346", size=4, parse = TRUE, hjust = 0) + 
  geom_text(data = d.eqn.80, aes(x = GL, y = BF, label = eqn),color="#6956E6", size=4, parse = TRUE, hjust = 0) + 
  geom_text(data = d.eqn.60, aes(x = GL, y = BF, label = eqn),color="#323189", size=4, parse = TRUE, hjust = 0) +
  labs(x="dGLS",y="ln(BF)") 

scat


ggsave(paste("Alldata",h,"scatter_lrmb.pdf",sep="_"), plot=scat,width = 8, height = 6, units = "in", device = 'pdf',bg = "transparent")




# Linear Scatter TvS -----------------------------------------------------------------------------------------------------------------------------------------------------

dataset <- "All"

tt <- paste(dataset,"Toxicofera v Scleroglossa",sep = " - ")


#magma(5)
# https://www.hexcolortool.com/#b5367e
quartz()

ds <- "Title"
h <- 'TvS'
cc <- c(  "midnightblue","slateblue1", "#B63679FF", "#FB8861FF", "#FCFDBFFF")
colorby <- "Subsample"
# Set graph stuff 


# dataset 
df <- df.m %>% filter(variable == h)

x.tic <- seq(-80,80,20) 
y.tic <- seq(-150,150,50)

# position of equation on chart
#xe <- mean(range(df$GL))
xe <- -40

d.eqn <- data.frame(GL = xe, BF = 160, eqn = lm_eqn(df$GL,df$BF))
d.eqn.80 <- data.frame(GL = xe, BF = 150, eqn = lm_eqn(df$GL80m,df$BF80m))
d.eqn.60 <- data.frame(GL = xe, BF = 140, eqn = lm_eqn(df$GL60m,df$BF60m))

scat <- ggplot(df) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80m, BF80m),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60m, BF60m),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.5,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80m,BF80m), method=lm,formula = y ~ x,size=0.5,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60m,BF60m), method=lm,formula = y ~ x,size=0.5,se=F, color="#323189",fullrange = T) +
  mytheme +
  scale_y_continuous(breaks = y.tic) + 
  scale_x_continuous(breaks = x.tic)+
  labs(x="dGLS",y="ln(BF)") + 
  geom_text(data = d.eqn, aes(x = GL, y = BF, label = eqn),color="#830346", size=2, parse = TRUE, hjust = 0) + 
  geom_text(data = d.eqn.80, aes(x = GL, y = BF, label = eqn),color="#6956E6", size=2, parse = TRUE, hjust = 0) + 
  geom_text(data = d.eqn.60, aes(x = GL, y = BF, label = eqn),color="#323189", size=2, parse = TRUE, hjust = 0) +
  labs(x="dGLS",y="ln(BF)") 

scat


ggsave(paste("Alldata",h,"scatter_lrm.pdf",sep="_"), plot=scat,width = 4, height = 3, units = "in", device = 'pdf',bg = "transparent")



# TVS mb -----------------------------------------------------------------------------------------------------------------------------------------------------

# dataset 
df <- df.mb %>% filter(variable == h)

x.tic <- seq(-80,80,20) 
y.tic <- seq(-150,150,50)

# position of equation on chart
#xe <- mean(range(df$GL))
xe <- -40

d.eqn <- data.frame(GL = xe, BF = 160, eqn = lm_eqn(df$GL,df$BF))
d.eqn.80 <- data.frame(GL = xe, BF = 150, eqn = lm_eqn(df$GL80mb,df$BF80mb))
d.eqn.60 <- data.frame(GL = xe, BF = 140, eqn = lm_eqn(df$GL60mb,df$BF60mb))

scat <- ggplot(df) + 
  geom_point(aes(GL,BF),alpha=0.5,size=0.5,color="#B63679FF") + 
  geom_point(aes(GL80mb, BF80mb),alpha=0.5,size=0.5,color="#826fff") + 
  geom_point(aes(GL60mb, BF60mb),alpha=0.5,size=0.5,color="#4B4AA2") + 
  geom_smooth(aes(GL,BF), method=lm, formula = y ~ x,size=0.5,se=F, color="#830346",fullrange = T) +
  geom_smooth(aes(GL80mb,BF80mb), method=lm,formula = y ~ x,size=0.5,se=F,color="#6956E6",fullrange = T) +
  geom_smooth(aes(GL60mb,BF60mb), method=lm,formula = y ~ x,size=0.5,se=F, color="#323189",fullrange = T) +
  mytheme +
  scale_y_continuous(breaks = y.tic) + 
  scale_x_continuous(breaks = x.tic)+
  labs(x="dGLS",y="ln(BF)") + 
  geom_text(data = d.eqn, aes(x = GL, y = BF, label = eqn),color="#830346", size=4, parse = TRUE, hjust = 0) + 
  geom_text(data = d.eqn.80, aes(x = GL, y = BF, label = eqn),color="#6956E6", size=4, parse = TRUE, hjust = 0) + 
  geom_text(data = d.eqn.60, aes(x = GL, y = BF, label = eqn),color="#323189", size=4, parse = TRUE, hjust = 0) +
  labs(x="dGLS",y="ln(BF)") 

scat


ggsave(paste("Alldata",h,"scatter_lrmb.pdf",sep="_"), plot=scat,width = 8, height = 6, units = "in", device = 'pdf',bg = "transparent")

## Depreciated 


# Mash data - another way
# Add column of Sample 
# all['Sample'] = 'all'
# all.80m['Sample'] = '80m'
# all.80mb['Sample'] = '80mb'
# all.60m['Sample'] = '60m'
# all.60mb['Sample'] = '60mb'
# # Mash
# x <- full_join(all.60m,all.80m, by=c('Locus', "variable", "GL", "BF", "dataSet", "diff"))
# y <- full_join(x,all, by=c('Locus', "variable", "GL", "BF", "dataSet", "diff"))
# z <- y %>% mutate(Subsample = coalesce(Sample.x, Sample.y,Sample))%>% select(-c(Sample.y, Sample.x, Sample))
# df.og <- z %>% filter(variable == h)
# scat.og <- ggplot(df.og, aes(x=GL,y=BF)) + 
#   geom_point(alpha=0.75, aes(color=Subsample), size=.75) + 
#   theme_bw() + theme(panel.border = element_blank(), legend.position = "none") + 
#   scale_color_manual(values=c("midnightblue","slateblue1", "#B63679FF")) + 
#   ggtitle("og")
# quartz()
# scat.og

#rm(list=setdiff(ls(), c("z")))

# lm_eqn = function(xvals,yvals){
#   m = lm(yvals ~ xvals) #y ~ x
#   eq <- substitute(italic(y) == a + b %.% italic(x)*",  "~~italic(r)^2~"="~r2, 
#                    list(a = format(unname(coef(m)[1]), digits = 2), 
#                         b = format(unname(coef(m)[2]), digits = 2), 
#                         r2 = format(summary(m)$r.squared, digits = 2)))
#   as.character(as.expression(eq));
# }
# 

