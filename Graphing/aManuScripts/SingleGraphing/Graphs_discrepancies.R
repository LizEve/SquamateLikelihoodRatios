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

dataset <- "Streicher"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))


# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS")) %>%
  mutate_if(is.numeric,round,2)
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF")) %>%
  mutate_if(is.numeric,round,2)

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10
mLBF[7:10] <- mLBF[7:10]/2

# Transform datasets - stack two comparisons, and add column for bf and gl 
# Each dataset one hypothesis 

B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))


G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType")) 

# Number of loci
loci <- length(mLGL$Locus)

# Smoosh together again
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
a <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% 
  select(-c("supportType.b","supportType.g")) 

d <- 2
d <- 0.5
# grouping when both are 0
# when 1 is zero, counting as "same"
s2 <- a %>% dplyr::rename(Hypothesis=variable) %>% 
  mutate(sig2 = case_when(BF==0 & GL==0 ~ 'both.zero.zero',
                          between(BF,0,5) & between(GL,0,d) ~ 'both.neu.same',
                       between(BF,-5,0) & between(GL,-d,0) ~ 'both.neu.same',
                       between(BF,0,5) & between(GL,-d,0) ~ 'both.neu.opposite',
                       between(BF,-5,0) & between(GL,0,d) ~ 'both.neu.opposite',
                       BF >= 5 & GL >= d ~ 'both.sig.same',
                       BF <= -5 & GL <= -d ~ 'both.sig.same',
                       BF <= -5 & GL >= d ~ 'both.sig.opposite',
                       BF >= 5 & GL <= -d ~ 'both.sig.opposite',
                       BF <= -5 & between(GL,-d,0) ~ 'neuG.sigB.same',
                       BF >= 5 & between(GL,0,d) ~ 'neuG.sigB.same',
                       GL <= -d & between(BF,-5,0) ~ 'neuB.sigG.same',
                       GL >= d & between(BF,0,5) ~ 'neuB.sigG.same',
                       BF <= -5 & between(GL,0,d) ~ 'neuG.sigB.opposite',
                       BF >= 5 & between(GL,-d,0) ~ 'neuG.sigB.opposite',
                       GL <= -d & between(BF,0,5) ~ 'neuB.sigG.opposite',
                       GL >= d & between(BF,-5,0) ~ 'neuB.sigG.opposite')) %>% 
  separate(sig2,c("a","b","c"),remove=FALSE) %>% 
  mutate(d=case_when(c=='same' & b == 'sigB' ~ 'same, one sig',
                     c=='same' & b == 'sigG' ~ 'same, one sig',
                     c=='same' & b == 'neu' ~ 'same, both neu',
                     c=='same' & b == 'sig' ~ 'same, both sig',
                     c=='opposite' ~ 'opposite',
                     c=='zero' ~ 'zero')) %>%
  mutate(Hypothesis=factor(Hypothesis, levels=c("AIvSA","AIvSI","SAvSI","TvS","ALL"))) %>%
  mutate(c=factor(c,levels = c("zero","opposite","same"))) %>%
  mutate(h=case_when(Hypothesis=='TvS' ~ 'TvS',
                     Hypothesis!='TvS' ~ 'Tox'))


# Count loci that have BF and dGLS that are in the same direction of support and those in opposite, and those that are zero
same.diff2 <- s2 %>% dplyr::group_by(c,Hypothesis) %>% dplyr::summarise(n.loci=n(),.groups = 'keep') %>% 
  mutate(percent.loci= round((n.loci/loci)*100,digits = 2)) %>% 
 ungroup() %>% group_by(Hypothesis) %>%  mutate(a=sum(n.loci),b=sum(percent.loci)) 



# Same as above but add int "same" when sig conflicts
same.diff.more2 <- s2 %>% dplyr::group_by(h,d) %>% dplyr::summarise(n=n(),.groups = 'keep') %>% 
  ungroup() %>% group_by(h) %>%
  dplyr::mutate(percent.loci= round((n/sum(n))*100,digits = 2)) %>% 
  ungroup() %>% group_by(h) %>%  mutate(a=sum(n),b=sum(percent.loci)) %>%
  separate(d,c("a","b","c"),remove=FALSE) 

# Tally all loci where significance thresholds conflict - both sig and opposite. one sig and opp, one sig and same
sig.diff2 <- s2  %>% 
  dplyr::group_by(sig2,h) %>% dplyr::summarise(n=n(),.groups = 'keep') %>% 
  ungroup() %>% group_by(h) %>%
  dplyr::mutate(percent.loci= round((n/sum(n))*100,digits = 2)) %>% 
  ungroup() %>% group_by(h) %>%  mutate(persum=sum(percent.loci)) %>%
  separate(sig2,c("a","b","c"),remove=FALSE) %>% 
  filter(!sig2 %in% c('both.neu.same','both.neu.opposite','both.zero.zero','both.sig.same'))


#save(overall,parsed,overall2,parsed2,z,z2, file=paste("Calcs_",dataset,".RData",sep=""))
# Graph Function-----------------------------------------------------------------------------------------------------------------------------------------------------
#position=position_dodge(), 
BarPlot <- function(df,H,yval,f,cc,tt){
  gg <- ggplot(df, aes(x=H,y=yval)) + 
    geom_bar(stat='identity', color="black",aes(fill=f),position=position_dodge()) +
    scale_fill_manual(values=cc, drop=FALSE)+
    theme_classic() + 
    theme(
      axis.text = element_text(size=12, color="black"),
      text = element_text(size=14),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5, vjust=0, size=14),
      legend.background = element_rect(colour = "transparent",fill = "transparent"),
      legend.title = element_text(size=12)) +
    labs(x='Direction of support',y='Percent of loci',fill="") +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    coord_cartesian(ylim=c(0,100)) + scale_y_continuous(breaks = seq(0,100,10)) +
    ggtitle(tt) 
  #+ scale_fill_grey(start=0.1,end = 0.9)
  return(gg)
}
setwd("/Users/ChatNoir/Projects/Squam/Graphs/ThresholdGraphs")
quartz()
cc <- viridis(4,direction = -1)
fig <- BarPlot(sig.diff2,sig.diff2$c,sig.diff2$percent.loci,sig.diff2$b,cc,paste("Conflict in support strength",dataset,sep=" ")) +
  coord_cartesian(ylim=c(0,35)) + scale_y_continuous(breaks = seq(0,35,5)) + 
  facet_wrap( ~ sig.diff2$h, strip.position = c("bottom")) + 
  scale_fill_manual(values=cc, drop=FALSE, labels=c("Both strong","BF strong, dGLS not", "dGLS strong, BF not"))
fig
# facet - , strip.position = c("bottom")
ggsave(paste(dataset,"discrepances.pdf",sep="_"), plot=fig,width = 6, height =4, units = "in", device = 'pdf',bg = "transparent")

cc <- viridis(3,direction = -1)
fig <- BarPlot(same.diff2,same.diff2$Hypothesis,same.diff2$percent.loci,same.diff2$c,
               cc,paste("Conflict in support strength",dataset,sep=" ")) + 
  labs(x='Hypothesis comparison',y='Percent of loci',fill="Direction \nof support")
fig
# facet - , strip.position = c("bottom")
ggsave(paste(dataset,"discrepances_general.pdf",sep="_"), plot=fig,width = 6, height =4, units = "in", device = 'pdf',bg = "transparent")

#------------------
# s1 <- a %>% dplyr::rename(Hypothesis=variable) %>% 
#   mutate(sig0.5 = case_when(between(BF,0,5) & between(GL,0,0.5) ~ 'neu.same',
#                             between(BF,-5,0) & between(GL,-0.5,0) ~ 'neu.same',
#                             between(BF,0,5) & between(GL,-0.5,0) ~ 'neu.opp',
#                             between(BF,-5,0) & between(GL,0,0.5) ~ 'neu.opp',
#                             BF >= 5 & GL >= 0.5 ~ 'sig.same',
#                             BF <= -5 & GL <= -0.5 ~ 'sig.same',
#                             BF <= -5 & GL >= 0.5 ~ 'sig.opp',
#                             BF >= 5 & GL <= -0.5 ~ 'sig.opp',
#                             BF <= -5 & between(GL,-0.5,0) ~ 'neuG.sigB.same',
#                             BF >= 5 & between(GL,0,0.5) ~ 'neuG.sigB.same',
#                             GL <= -0.5 & between(BF,-5,0) ~ 'neuB.sigG.same',
#                             GL >= 0.5 & between(BF,0,5) ~ 'neuB.sigG.same',
#                             BF <= -5 & between(GL,0,0.5) ~ 'neuG.sigB.opp',
#                             BF >= 5 & between(GL,-0.5,0) ~ 'neuG.sigB.opp',
#                             GL <= -0.5 & between(BF,0,5) ~ 'neuB.sigG.opp',
#                             GL >= 0.5 & between(BF,-5,0) ~ 'neuB.sigG.opp'))
# 
# s <- left_join(s1,s2,by=c("Locus","Hypothesis", "GL", "BF"))
# # Calculate metrics for all loci, all loci is divided by 4 times loci to account for all hypotheses
# x <- s %>% group_by(sig0.5) %>% dplyr::summarise(Hypothesis='ALL',n.loci=n()) %>%
#   mutate(percent.loci =  round((n.loci/(loci*4))*100,digits=2))
# y <- s %>% group_by(sig0.5,Hypothesis) %>% dplyr::summarise(n.loci=n()) %>% 
#   arrange(sig0.5,Hypothesis) %>% 
#   mutate(percent.loci= round((n.loci/loci)*100,digits = 2))
# # This will produce an error because it removes factors, this is OKAY.
# z <- bind_rows(x,y) %>% separate(sig0.5,c("a","b","c")) %>% 
#   mutate(Hypothesis=factor(Hypothesis, levels=c("AIvSA","AIvSI","SAvSI","TvS","ALL")))
# 




### This will produce an error because it removes factors, this is OKAY.
##z2 <- y2 %>% separate(sig2,c("a","b","c"),remove=FALSE) %>% 
##  mutate(Hypothesis=factor(Hypothesis, levels=c("AIvSA","AIvSI","SAvSI","TvS","ALL"))) 
##
##z2 <- y2 %>% mutate(Hypothesis=factor(Hypothesis, levels=c("AIvSA","AIvSI","SAvSI","TvS","ALL"))) 
##
##
##overall <- z %>% filter(is.na(c)) %>% select(-c) %>% unite("conflict",a:b,remove=FALSE) 
##parsed <- z %>% filter(!is.na(c)) %>% unite("conflict",c(a,b,c),remove=FALSE)
##
##overall2 <- z2 %>% filter(is.na(c)) %>% select(-c) %>% unite("conflict",a:b,remove=FALSE)
##parsed2 <- z2 %>% filter(!is.na(c)) %>% unite("conflict",c(a,b,c),remove=FALSE)
##



















