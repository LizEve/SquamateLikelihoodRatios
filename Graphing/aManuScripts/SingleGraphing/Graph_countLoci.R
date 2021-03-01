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
library(forcats)
library(DescTools)

load("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/Calcs_countLociALL.RData")
load("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/Calcs_countLociALL2.RData")

# set up to load both datasets to graph if you want. i haven't dont that yet. 

#--- Currently using this data structure ---------------------------------------------------
x <- c %>% mutate(Hsig = case_when(variable == "AIvSA" & value == "Pos Sig" ~ "AI",
                              variable == "AIvSA" & value == "Neg Sig" ~ "SA",
                              variable == "AIvSI" & value == "Pos Sig" ~ "AI",
                              variable == "AIvSI" & value == "Neg Sig" ~ "SI",
                              variable == "SAvSI" & value == "Pos Sig" ~ "SA",
                              variable == "SAvSI" & value == "Neg Sig" ~ "SI",
                              variable == "TvS" & value == "Pos Sig" ~ "TX",
                              variable == "TvS" & value == "Neg Sig" ~ "SC",
                              value == "Not Sig" ~ "Not Sig"))  %>% 
  mutate(Hpos = case_when(variable == "AIvSA" & value == "Pos" ~ "AI",
                          variable == "AIvSA" & value == "Neg" ~ "SA",
                          variable == "AIvSI" & value == "Pos" ~ "AI",
                          variable == "AIvSI" & value == "Neg" ~ "SI",
                          variable == "SAvSI" & value == "Pos" ~ "SA",
                          variable == "SAvSI" & value == "Neg" ~ "SI",
                          variable == "TvS" & value == "Pos" ~ "TX",
                          variable == "TvS" & value == "Neg" ~ "SC",
                          value=="Zero" ~ "Zero")) %>%
  mutate(Hpos=factor(Hpos, levels=c("AI","SA","SI","TX","SC","Not Sig","Zero")))  %>%
  mutate(Hsig=factor(Hsig, levels=c("AI","SA","SI","TX","SC","Not Sig","Zero"))) %>% 
  mutate(variable=factor(variable, levels=c("AIvSA","AIvSI","SAvSI","TvS")))
g <- x %>% filter(metric=="GL") %>% dplyr::rename(GL.percent.loci=Percent.Loci) %>% select(-c("metric","info"))
b <- x %>% filter(metric=="BF") %>% dplyr::rename(BF.percent.loci=Percent.Loci) %>% select(-c("metric","info"))
all <- inner_join(g,b,by=c("Data.Set","value","variable","Hpos","Hsig"), suffix = c(".g", ".b"))
all.sig <- all %>% filter(!is.na(Hsig)) %>% select(-c("Hpos")) 
all.pos <- all %>% filter(!is.na(Hpos)) %>% select(-c("Hsig")) 


x2 <- c2 %>% mutate(Hsig = case_when(variable == "AIvSA" & value == "Pos Sig" ~ "AI",
                                   variable == "AIvSA" & value == "Neg Sig" ~ "SA",
                                   variable == "AIvSI" & value == "Pos Sig" ~ "AI",
                                   variable == "AIvSI" & value == "Neg Sig" ~ "SI",
                                   variable == "SAvSI" & value == "Pos Sig" ~ "SA",
                                   variable == "SAvSI" & value == "Neg Sig" ~ "SI",
                                   variable == "TvS" & value == "Pos Sig" ~ "TX",
                                   variable == "TvS" & value == "Neg Sig" ~ "SC",
                                   value == "Not Sig" ~ "Not Sig"))  %>% 
  mutate(Hpos = case_when(variable == "AIvSA" & value == "Pos" ~ "AI",
                          variable == "AIvSA" & value == "Neg" ~ "SA",
                          variable == "AIvSI" & value == "Pos" ~ "AI",
                          variable == "AIvSI" & value == "Neg" ~ "SI",
                          variable == "SAvSI" & value == "Pos" ~ "SA",
                          variable == "SAvSI" & value == "Neg" ~ "SI",
                          variable == "TvS" & value == "Pos" ~ "TX",
                          variable == "TvS" & value == "Neg" ~ "SC",
                          value=="Zero" ~ "Zero")) %>%
  mutate(Hpos=factor(Hpos, levels=c("AI","SA","SI","TX","SC","Not Sig","Zero")))  %>%
  mutate(Hsig=factor(Hsig, levels=c("AI","SA","SI","TX","SC","Not Sig","Zero"))) %>% 
  mutate(variable=factor(variable, levels=c("AIvSA","AIvSI","SAvSI","TvS")))
g2 <- x2 %>% filter(metric=="GL") %>% dplyr::rename(GL.percent.loci=Percent.Loci) %>% select(-c("metric","info"))
b2 <- x2 %>% filter(metric=="BF") %>% dplyr::rename(BF.percent.loci=Percent.Loci) %>% select(-c("metric","info"))
all2 <- inner_join(g2,b2,by=c("Data.Set","value","variable","Hpos","Hsig"), suffix = c(".g", ".b"))
all.sig2 <- all2 %>% filter(!is.na(Hsig)) %>% select(-c("Hpos")) 
all.pos2 <- all2 %>% filter(!is.na(Hpos)) %>% select(-c("Hsig")) 


#--- Graph Function---------------------------------------------------

BarPlot <- function(df,H,yval,f,cc,tt){
  gg <- ggplot(df, aes(x=H,y=yval,fill=f)) + 
    geom_bar(stat='identity', position=position_dodge(),color="black") +
    scale_color_manual(values=cc, drop=FALSE) + scale_fill_manual(values=cc, drop=FALSE) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=12, color="black"),
      text = element_text(size=14),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.3, vjust=-4, size=14),
      legend.background = element_rect(colour = "transparent",fill = "transparent"),
      legend.title = element_text(size=12)) +
    #labs(x="Hypotheses compared",y='Percent of loci',fill="Hypothesis \n supported") +
    labs(x="",y='',fill="Hypothesis \n supported") +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    coord_cartesian(ylim=c(0,100)) + scale_y_continuous(breaks = seq(0,100,20)) +
    ggtitle(tt) + 
    geom_text(label = paste(round(yval,0),"%",sep=''),
                            position = position_dodge(0.9),
                            vjust = -0.5,size=3)
  #+ scale_fill_grey(start=0.1,end = 0.9,drop=FALSE)
  return(gg)
}


#--- Graph Data---------------------------------------------------

setwd("/Users/ChatNoir/Projects/Squam/Graphs/PhyloSupportGraphs")

#preatty coooolars
cc <- c(viridis(5), "gray90", "white")
#cc <- viridis(7, option="C")
#cc <- ColToGrey(cc)
# Subset data 
dataset <- "Burbrink"
ap <- all.pos %>% filter(Data.Set==dataset)
as <- all.sig %>% filter(Data.Set==dataset)

# BF

bp <- BarPlot(ap,ap$variable,ap$BF.percent.loci,ap$Hpos,cc,"Positive ln(BF) support")
bs <- BarPlot(as,as$variable,as$BF.percent.loci,as$Hsig,cc,"Strong ln(BF) support") 

# GL
gp <- BarPlot(ap,ap$variable,ap$GL.percent.loci,ap$Hpos,cc,"Positive dGLS support")
gs <- BarPlot(as,as$variable,as$GL.percent.loci,as$Hsig,cc,"Strong dGLS support") 

p <- ggarrange(bp,bs,gp,gs, ncol=2, nrow=2,align="v", common.legend = TRUE, legend = "right")

fig <- annotate_figure(p,
                bottom = text_grob("Hypotheses compared", size = 16),
                left = text_grob("Percent of loci",  size = 16, rot = 90))
#fig.lab = "Figure 1", fig.lab.face = "bold", fig.lab.size = 18
fig

ggsave(paste(dataset,"bar_support.pdf",sep="_"), plot=fig,width = 12, height = 10, units = "in", device = 'pdf',bg = "transparent")

#--- dGLS = 2 ---------------------------------------------------


ap <- all.pos2 %>% filter(Data.Set==dataset)
as <- all.sig2 %>% filter(Data.Set==dataset)

# BF
bp <- BarPlot(ap,ap$variable,ap$BF.percent.loci,ap$Hpos,cc,"Positive ln(BF) support")
bs <- BarPlot(as,as$variable,as$BF.percent.loci,as$Hsig,cc,"Strong ln(BF) support") 

# GL
gp <- BarPlot(ap,ap$variable,ap$GL.percent.loci,ap$Hpos,cc,"Positive dGLS support")
gs <- BarPlot(as,as$variable,as$GL.percent.loci,as$Hsig,cc,"Strong dGLS support") 

p <- ggarrange(bp,bs,gp,gs, ncol=2, nrow=2,align="v", common.legend = TRUE, legend = "right")

fig <- annotate_figure(p,
                       bottom = text_grob("Hypotheses compared", size = 16),
                       left = text_grob("Percent of loci",  size = 16, rot = 90))
#fig.lab = "Figure 1", fig.lab.face = "bold", fig.lab.size = 18
fig


#ggsave(paste(dataset,"bar_support2.pdf",sep="_"), plot=fig,width = 12, height = 10, units = "in", device = 'pdf',bg = "transparent")

#--- Different data structure ---------------------------------------------------
#y <- c %>% mutate(Hypoth.Supported = case_when(variable == "AIvSA" & value == "Pos Sig" ~ "AI Sig",
#                                               variable == "AIvSA" & value == "Neg Sig" ~ "SA Sig",
#                                               variable == "AIvSI" & value == "Pos Sig" ~ "AI Sig",
#                                               variable == "AIvSI" & value == "Neg Sig" ~ "SI Sig",
#                                               variable == "SAvSI" & value == "Pos Sig" ~ "SA Sig",
#                                               variable == "SAvSI" & value == "Neg Sig" ~ "SI Sig",
#                                               variable == "TvS" & value == "Pos Sig" ~ "TX Sig",
#                                               variable == "TvS" & value == "Neg Sig" ~ "SC Sig",
#                                               variable == "AIvSA" & value == "Pos" ~ "AI Pos",
#                                               variable == "AIvSA" & value == "Neg" ~ "SA Pos",
#                                               variable == "AIvSI" & value == "Pos" ~ "AI Pos",
#                                               variable == "AIvSI" & value == "Neg" ~ "SI Pos",
#                                               variable == "SAvSI" & value == "Pos" ~ "SA Pos",
#                                               variable == "SAvSI" & value == "Neg" ~ "SI Pos",
#                                               variable == "TvS" & value == "Pos" ~ "TX Pos",
#                                               variable == "TvS" & value == "Neg" ~ "SC Pos", 
#                                               value == "Zero" ~ "Zero",
#                                               value == "Not Sig" ~ "Not Sig")) %>%
#  mutate(Hypoth.Supported=factor(Hypoth.Supported, levels=c("AI Sig","AI Pos",
#                                                            "SA Sig","SA Pos",
#                                                            "SI Sig","SI Pos",
#                                                            "TX Sig","TX Pos",
#                                                            "SC Sig","SC Pos",
#                                                            "Not Sig","Zero")))
#
#
#g <- y %>% filter(metric=="GL") %>% dplyr::rename(GL.percent.loci=Percent.Loci) %>% select(-c("metric","info"))
#b <- y %>% filter(metric=="BF") %>% dplyr::rename(BF.percent.loci=Percent.Loci) %>% select(-c("metric","info"))
#all <- inner_join(g,b,by=c("Data.Set","value","variable","Hypoth.Supported"), suffix = c(".g", ".b"))
#

