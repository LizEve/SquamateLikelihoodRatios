rm(list=ls())
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)

library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library("viridis")
library(forcats)
library(DescTools)

load("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/Calcs_countLociALL.RData")
load("/Users/ChatNoir/Projects/Squam/scripts_ch1/Graphing/DataFiles/Calcs_countLociALL2.RData")

# set up to load both datasets to graph if you want. i haven't dont that yet. 


#--- Currently using this data structure - but not this data ---------------------------------------------------
#x <- c %>%  dplyr::mutate(Hsig = case_when(variable == "AIvSA" & value == "Pos Sig" ~ "AI",
#                              variable == "AIvSA" & value == "Neg Sig" ~ "SA",
#                              variable == "AIvSI" & value == "Pos Sig" ~ "AI",
#                              variable == "AIvSI" & value == "Neg Sig" ~ "SI",
#                              variable == "SAvSI" & value == "Pos Sig" ~ "SA",
#                              variable == "SAvSI" & value == "Neg Sig" ~ "SI",
#                              variable == "TvS" & value == "Pos Sig" ~ "TX",
#                              variable == "TvS" & value == "Neg Sig" ~ "SC",
#                              value == "Not Sig" ~ "Not Sig"))  %>% 
#  dplyr::mutate(Hpos = case_when(variable == "AIvSA" & value == "Pos" ~ "AI",
#                          variable == "AIvSA" & value == "Neg" ~ "SA",
#                          variable == "AIvSI" & value == "Pos" ~ "AI",
#                          variable == "AIvSI" & value == "Neg" ~ "SI",
#                          variable == "SAvSI" & value == "Pos" ~ "SA",
#                          variable == "SAvSI" & value == "Neg" ~ "SI",
#                          variable == "TvS" & value == "Pos" ~ "TX",
#                          variable == "TvS" & value == "Neg" ~ "SC",
#                          value=="Zero" ~ "Zero")) %>%
#  dplyr::mutate(Hpos=factor(Hpos, levels=c("AI","SA","SI","TX","SC","Not Sig","Zero")))  %>%
#  dplyr::mutate(Hsig=factor(Hsig, levels=c("AI","SA","SI","TX","SC","Not Sig","Zero"))) %>% 
#  dplyr::mutate(variable=factor(variable, levels=c("AIvSA","AIvSI","SAvSI","TvS")))
#
#
#g <- x %>%  dplyr::filter(metric=="GL") %>% dplyr::rename(GL.percent.loci=Percent.Loci) %>% select(-c("metric","info"))
#b <- x %>%  dplyr::filter(metric=="BF") %>% dplyr::rename(BF.percent.loci=Percent.Loci) %>% select(-c("metric","info"))
#all <- inner_join(g,b,by=c("Data.Set","value","variable","Hpos","Hsig"), suffix = c(".g", ".b"))
#all.sig <- all %>%  dplyr::filter(!is.na(Hsig)) %>% select(-c("Hpos")) 
#all.pos <- all %>%  dplyr::filter(!is.na(Hpos)) %>% select(-c("Hsig")) 

#--- Currently using this data structure 2---------------------------------------------------

x2 <- c2 %>% mutate(Hsig = case_when(variable == "AIvSA" & value == "Pos Sig" ~ "AI",
                                   variable == "AIvSA" & value == "Neg Sig" ~ "SA",
                                   variable == "AIvSI" & value == "Pos Sig" ~ "AI",
                                   variable == "AIvSI" & value == "Neg Sig" ~ "SI",
                                   variable == "SAvSI" & value == "Pos Sig" ~ "SA",
                                   variable == "SAvSI" & value == "Neg Sig" ~ "SI",
                                   variable == "TvS" & value == "Pos Sig" ~ "TX",
                                   variable == "TvS" & value == "Neg Sig" ~ "SC",
                                   value == "Not Sig" ~ "Weak or neutral"))  %>% 
  mutate(Hpos = case_when(variable == "AIvSA" & value == "Pos" ~ "AI",
                          variable == "AIvSA" & value == "Neg" ~ "SA",
                          variable == "AIvSI" & value == "Pos" ~ "AI",
                          variable == "AIvSI" & value == "Neg" ~ "SI",
                          variable == "SAvSI" & value == "Pos" ~ "SA",
                          variable == "SAvSI" & value == "Neg" ~ "SI",
                          variable == "TvS" & value == "Pos" ~ "TX",
                          variable == "TvS" & value == "Neg" ~ "SC",
                          value=="Zero" ~ "Zero")) %>%
  mutate(Hpos=factor(Hpos, levels=c("AI","SA","SI","TX","SC","Weak or neutral","Zero")))  %>%
  mutate(Hsig=factor(Hsig, levels=c("AI","SA","SI","TX","SC","Weak or neutral","Zero"))) %>% 
  mutate(variable=factor(variable, levels=c("AIvSA","AIvSI","SAvSI","TvS")))

pos <- x2 %>% filter(!is.na(Hpos)) %>% select(-c(Hsig))
sig <- x2 %>% filter(!is.na(Hsig)) %>% select(-c(Hpos))


TS.pos <- pos %>% filter(variable=="TvS") %>% 
  droplevels() %>% 
  mutate(Hpos=factor(Hpos, levels=c("TX","SC","Zero","Weak or neutral"))) %>% 
  complete(Hpos,Data.Set,variable,metric) %>% 
  replace_na(list(Percent.Loci = 0)) %>%
  filter(Hpos!="Weak or neutral")

G.TS.pos <- TS.pos %>% filter(metric=="GL")
B.TS.pos <- TS.pos %>% filter(metric=="BF")

TS.sig <- sig %>% filter(variable=="TvS") %>% 
  droplevels() %>% 
  mutate(Hsig=factor(Hsig, levels=c("TX","SC","Zero","Weak or neutral"))) %>% 
  complete(Hsig,Data.Set,variable) %>% 
  replace_na(list(Percent.Loci = 0))

G.TS.sig <- TS.sig %>% filter(metric=="GL")
B.TS.sig <- TS.sig %>% filter(metric=="BF")

TX.pos <- pos %>% filter(variable!="TvS") %>% 
  droplevels() %>%
  mutate(Hpos=factor(Hpos, levels=c("AI","SA","SI","Zero","Weak or neutral"))) %>% 
  complete(Hpos,Data.Set,variable,metric) %>% 
  replace_na(list(Percent.Loci = 0)) %>% 
  filter(!(Hpos=="AI"& variable=="SAvSI")) %>%
  filter(!(Hpos=="SA"& variable=="AIvSI")) %>%
  filter(!(Hpos=="SI"& variable=="AIvSA")) %>%
  filter(Hpos!="Weak or neutral")

G.TX.pos <- TX.pos %>% filter(metric=="GL")
B.TX.pos <- TX.pos %>% filter(metric=="BF")

TX.sig <- sig %>% filter(variable!="TvS") %>% 
  droplevels() %>%
  mutate(Hsig=factor(Hsig, levels=c("AI","SA","SI","Zero","Weak or neutral"))) %>% 
  complete(Hsig,Data.Set,variable,metric) %>% 
  replace_na(list(Percent.Loci = 0))%>% 
  filter(!(Hsig=="AI"& variable=="SAvSI")) %>%
  filter(!(Hsig=="SA"& variable=="AIvSI")) %>%
  filter(!(Hsig=="SI"& variable=="AIvSA")) 

G.TX.sig <- TX.sig %>% filter(metric=="GL")
B.TX.sig <- TX.sig %>% filter(metric=="BF")

#--- Graph Function---------------------------------------------------
#position=position_dodge(), 
BarPlot <- function(df,H,yval,f,cc,tt){
  gg <- ggplot(df, aes(x=H,y=yval,fill=f)) + 
    geom_bar(stat='identity', color="black", position=position_dodge()) +
    scale_color_manual(values=cc, drop=FALSE) + scale_fill_manual(values=cc, drop=FALSE) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=8, color="black"),
      axis.text.x = element_text(angle = 45, hjust=1),
      text = element_text(size=10),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5, vjust=0, size=10),
      legend.background = element_rect(colour = "transparent",fill = "transparent"),
      legend.title = element_text(size=10)) +
    #labs(x="Hypotheses compared",y='Percent of loci',fill="Hypothesis \n supported") +
    labs(x="",y='',fill="Hypothesis \n supported") +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    coord_cartesian(ylim=c(0,100)) + scale_y_continuous(breaks = seq(0,100,10)) +
    ggtitle(tt)
  #+ scale_fill_grey(start=0.1,end = 0.9,drop=FALSE)
  return(gg)
}

#--- Colors and Folder---------------------------------------------------

setwd("/Users/ChatNoir/Projects/Squam/Graphs/PhyloSupportGraphs")
#preatty coooolars
cc <- c(viridis(5), "gray90", "white")
#cc <- viridis(7, option="C")
#cc <- ColToGrey(cc)
quartz()

#--- TX - Graph Data---------------------------------------------------

G.TX.pos
B.TX.pos
G.TX.sig
B.TX.sig


cc <- c(viridis(4)[1:3], "gray90", "white")
bp <- BarPlot(B.TX.pos,B.TX.pos$Data.Set,B.TX.pos$Percent.Loci,B.TX.pos$Hpos,cc,"a. ln(BF) - Positive support") + 
  facet_wrap(~ B.TX.pos$variable,ncol=4,strip.position = c("bottom")) + 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
bp
bs <- BarPlot(B.TX.sig,B.TX.sig$Data.Set,B.TX.sig$Percent.Loci,B.TX.sig$Hsig,cc,"b. ln(BF) - Strong support") + 
  facet_wrap(~ B.TX.sig$variable,ncol=4,strip.position = c("bottom"))+ 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
bs
gp <- BarPlot(G.TX.pos,G.TX.pos$Data.Set,G.TX.pos$Percent.Loci,G.TX.pos$Hpos,cc,"c. dGLS - Positive support") + 
  facet_wrap(~ G.TX.pos$variable,ncol=4,strip.position = c("bottom"))+ 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
gs <- BarPlot(G.TX.sig,G.TX.sig$Data.Set,G.TX.sig$Percent.Loci,G.TX.sig$Hsig,cc,"d. dGLS - Strong support") + 
  facet_wrap(~ G.TX.sig$variable,ncol=4,strip.position = c("bottom"))+ 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))


p <- ggarrange(bp,bs,gp,gs, ncol=2, nrow=2,align="v", common.legend = TRUE, legend = "bottom")

fig <- annotate_figure(p,left = text_grob("Percent of loci",  size = 12, vjust = 0.5, rot = 90))

#fig <- annotate_figure(p,
#                bottom = text_grob("Hypotheses compared", size = 12,vjust = -0.5),
#                left = text_grob("Percent of loci",  size = 12, vjust = 0.5, rot = 90))
#fig.lab = "Figure 1", fig.lab.face = "bold", fig.lab.size = 18
fig

#ggsave("Fig6_All_bar_support_TX.pdf", plot=fig,width = 8, height = 8, units = "in", device = 'pdf',bg = "transparent")

#--- TvS ---------------------------------------------------
G.TS.pos
B.TS.pos
G.TS.sig
B.TS.sig

cc <- c(viridis(5)[4:5], "gray90", "white")
bp <- BarPlot(B.TS.pos,B.TS.pos$Data.Set,B.TS.pos$Percent.Loci,B.TS.pos$Hpos,cc,"a. ln(BF) - Positive support") + 
  facet_wrap(~ B.TS.pos$variable,ncol=4,strip.position = c("bottom"))
bs <- BarPlot(B.TS.sig,B.TS.sig$Data.Set,B.TS.sig$Percent.Loci,B.TS.sig$Hsig,cc,"b. ln(BF) - Strong support") + 
  facet_wrap(~ B.TS.sig$variable,ncol=4,strip.position = c("bottom"))

gp <- BarPlot(G.TS.pos,G.TS.pos$Data.Set,G.TS.pos$Percent.Loci,G.TS.pos$Hpos,cc,"c. dGLS - Positive support") + 
  facet_wrap(~ G.TS.pos$variable,ncol=4,strip.position = c("bottom"))
gs <- BarPlot(G.TS.sig,G.TS.sig$Data.Set,G.TS.sig$Percent.Loci,G.TS.sig$Hsig,cc,"d. dGLS - Strong support") + 
  facet_wrap(~ G.TS.sig$variable,ncol=4,strip.position = c("bottom"))


p <- ggarrange(bp,bs,gp,gs, ncol=2, nrow=2,align="v", common.legend = TRUE, legend = "bottom")

fig <- annotate_figure(p,left = text_grob("Percent of loci",  size = 12, vjust = 0.5, rot = 90))

#fig.lab = "Figure 1", fig.lab.face = "bold", fig.lab.size = 18
fig

ggsave("SIFigx_All_bar_support_TS.pdf", plot=fig,width = 6, height = 6, units = "in", device = 'pdf',bg = "transparent")




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

