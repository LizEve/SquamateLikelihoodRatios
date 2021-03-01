#rm(list=ls())

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
#rm(a,all.e,B,df,fig,G,mLBF,mLGL,s,support.direction,t.hold,t.hold.details)
dataset <- "Singhal"

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
#d <- 0.5
# grouping when both are 0
# when 1 is zero, counting as "same"
s <- a %>% dplyr::rename(Hypothesis=variable) %>% 
  mutate(sig = case_when(BF==0 & GL==0 ~ 'both.zero.zero',
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
  separate(sig,c("a","b","c"),remove=FALSE) %>% 
  mutate(d=case_when(c=='same' & b == 'sigB' ~ 'same, one sig',
                     c=='same' & b == 'sigG' ~ 'same, one sig',
                     c=='same' & b == 'neu' ~ 'same, both neu',
                     c=='same' & b == 'sig' ~ 'same, both sig',
                     c=='opposite' ~ 'opposite',
                     c=='zero' ~ 'zero')) %>%
  mutate(Hypothesis=factor(Hypothesis, levels=c("AIvSA","AIvSI","SAvSI","TvS"))) %>%
  mutate(c=factor(c,levels = c("zero","opposite","same"))) %>%
  mutate(h=case_when(Hypothesis=='TvS' ~ 'TvS',
                     Hypothesis!='TvS' ~ 'Tox')) %>% 
  mutate(thresh=case_when(a == 'neuB' ~ 'conflict',
                          a == 'neuG' ~ 'conflict',
                          sig == 'both.sig.opposite' ~ 'conflict',
                          sig == 'both.neu.opposite' ~ 'agree',
                          sig == 'both.neu.same' ~ 'agree',
                          sig == 'both.sig.same' ~ 'agree',
                          sig == 'both.zero.zero' ~ 'agree'))


# Count loci that have BF and dGLS that are in the same direction of support and those in opposite, and those that are zero
support.direction <- s %>% dplyr::group_by(c,Hypothesis,.drop=FALSE) %>% dplyr::summarise(n.loci=n(),.groups = 'keep') %>% 
  mutate(percent.loci= round((n.loci/loci)*100,digits = 2)) %>% 
 ungroup() %>% group_by(Hypothesis) %>%  mutate(a=sum(n.loci),b=sum(percent.loci)) 


# Count loci that disagree in terms of thresholds 
t.hold <- s %>% dplyr::group_by(thresh,Hypothesis,.drop=FALSE) %>% dplyr::summarise(n.loci=n(),.groups = 'keep') %>% 
  mutate(percent.loci= round((n.loci/loci)*100,digits = 2)) %>% 
  ungroup() %>% group_by(Hypothesis) %>%  mutate(a=sum(n.loci),b=sum(percent.loci)) %>%
  mutate(thresh=factor(thresh,levels = c("conflict","agree")))




# Tally all loci where significance thresholds conflict - both sig and opposite. one sig and opp, one sig and same
t.hold.details <- s  %>% filter(!sig %in% c('both.neu.same','both.neu.opposite','both.zero.zero','both.sig.same')) %>%
  separate(sig,c("a","b","c"),remove=FALSE) %>%
  unite("ab",a:b, remove=FALSE) %>%
  dplyr::group_by(ab,Hypothesis,.drop=FALSE) %>% dplyr::summarise(n=n(),.groups = 'keep') %>% 
  ungroup() %>% group_by(Hypothesis) %>%
  dplyr::mutate(percent.loci= round((n/loci)*100,digits = 2)) %>% 
  droplevels() %>% mutate(thresh=factor(ab,levels = c("both_sig","neuG_sigB","neuB_sigG")))


#save(overall,parsed,overall2,parsed2,z,z2, file=paste("Calcs_",dataset,".RData",sep=""))
# Graph Function-----------------------------------------------------------------------------------------------------------------------------------------------------
#position=position_dodge(), 
BarPlot <- function(df,H,yval,f,cc,tt,lll,tf){
  gg <- ggplot(df, aes(x=H,y=yval)) + 
    geom_bar(stat='identity', color="black",aes(fill=f),position=position_dodge()) +
    scale_fill_manual(values=cc, drop=FALSE,labels=lll,guide = guide_legend(reverse=tf))+
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
#quartz()

df <- support.direction
cc <- viridis(3,direction = -1)
fig <- BarPlot(df,df$Hypothesis,df$percent.loci,df$c,
               cc,paste("Conflicting ln(BF) and dGLS values -",dataset,sep=" "),
               c("zero", "conflict", "agree"), TRUE) + 
  labs(x='Hypothesis comparison',y='Percent of loci',fill="Direction \nof support") + 
  geom_text(aes(x=df$Hypothesis, y=df$percent.loci, 
                label = paste(round(df$percent.loci,0),"%",sep=''),
                group=df$c),
            position = position_dodge(width=1),
            vjust = -0.5,size=4)
#fig
# facet - , strip.position = c("bottom")
#ggsave(paste(dataset,"discrepances_general.pdf",sep="_"), plot=fig1,width = 6, height =4, units = "in", device = 'pdf',bg = "transparent")

df <- t.hold
cc <- viridis(3,direction = -1)
cc <- c("#21908CFF", "#440154FF")
fig <- BarPlot(df,df$Hypothesis,df$percent.loci,df$thresh,
               cc,paste("Strong support threshold conflicts -",dataset,sep=" "),
               c( "conflict","agree"),TRUE) + 
  labs(x='Hypothesis comparison',y='Percent of loci',fill="Direction \nof support") + 
  geom_text(aes(x=df$Hypothesis, y=df$percent.loci, 
                label = paste(round(df$percent.loci,0),"%",sep=''),
                group=df$thresh),
            position = position_dodge(width=1),
            vjust = -0.5,size=4)
#fig
#ggsave(paste(dataset,"discrepances_threshold.pdf",sep="_"), plot=fig,width = 6, height =4, units = "in", device = 'pdf',bg = "transparent")


df <- t.hold.details %>% ungroup() %>% select(c("Hypothesis","percent.loci","thresh"))
cc <- viridis(6,direction = -1)
fig <- BarPlot(df,df$Hypothesis,df$percent.loci,df$thresh,cc,paste("Conflicting strong support by metric",dataset,sep=" "),
               c("Both strong","BF strong", "dGLS strong"),TRUE) +
  coord_cartesian(ylim=c(0,55)) + scale_y_continuous(breaks = seq(0,55,5)) + 
  labs(x='Hypothesis comparison',y='Percent of loci',fill="Conflicts") + 
  geom_text(aes(x=Hypothesis, y=percent.loci, 
                label = paste(round(percent.loci,0),"%",sep=''),
                group=thresh),
            position = position_dodge(width=1),
            vjust = -0.5,size=4)
#fig
#ggsave(paste(dataset,"discrepances_threshold_conflict.pdf",sep="_"), plot=fig,width = 6, height =4, units = "in", device = 'pdf',bg = "transparent")


