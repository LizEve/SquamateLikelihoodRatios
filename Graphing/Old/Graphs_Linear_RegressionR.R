rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
library(praise)
praise()

# set dataset 

dataset <- "Reeder"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS")) %>%
  mutate_if(is.numeric,round,2)
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF")) %>%
  mutate_if(is.numeric,round,2)

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10
mLBF[7:10] <- mLBF[7:10]/2

B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))


G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))

# Smoosh together again
a <- bind_rows(B,G) 
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
all.exp <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% select(-c("supportType.b","supportType.g"))
all <- all.exp
# Disagree Data ------------------------------------------------------------------------------------------------------------------------------------------------------

all.d <- all %>% 
  mutate(diff = case_when(between(GL,-2,2) & BF > 5 ~ "BF strong dGLS neutral",
                          between(GL,-2,2) & BF < -5 ~ "BF strong dGLS neutral",
                          between(BF,-5,5) & GL > 2 ~ "dGLS strong BF neutral",
                          between(BF,-5,5) & GL < -2 ~ "dGLS strong BF neutral",
                          BF < -5 & GL > 2 ~ "opposite strong support",
                          BF > 5 & GL < -2 ~ "opposite strong support",
                          BF > 5 & GL > 2 ~ "agree",
                          BF < -5 & GL < -2 ~ "agree",
                          between(BF,-5,5) & between(GL,-2,2) ~ "agree"))

# Subset of Data ------------------------------------------------------------------------------------------------------------------------------------------------------


all.2x <- all.d %>% group_by(variable) %>% 
  filter(between(BF,-15,15) | 
           between(GL,-6,6)) # "significance levels (5,2) 5 + 5x 2
all.1x <- all.d %>% group_by(variable) %>% 
  filter(between(BF,-10,10) | 
           between(GL,-4,4)) # "significance levels (5,2) 5 + 5x 1
all.0x <- all.d %>% group_by(variable) %>% 
  filter(between(BF,-5,5) | 
           between(GL,-2,2)) # "significance levels (5,2) 5 + 5x 1


# get percentiles 

all.50l <- all.d %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.d$GL,0),quantile(all.d$GL, 0.5)) |
           between(BF,quantile(all.d$BF,0),quantile(all.d$BF, 0.5)))

all.50h <- all.d %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.d$GL,0.5),quantile(all.d$GL, 1)) |
           between(BF,quantile(all.d$BF,0.5),quantile(all.d$BF, 1)))

all.50m <- all.d %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.d$GL,0.25),quantile(all.d$GL, 0.75)) |
           between(BF,quantile(all.d$BF,0.25),quantile(all.d$BF, 0.75)))

all.60m <- all.d %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.d$GL,0.2),quantile(all.d$GL, 0.8)) | 
           between(BF,quantile(all.d$BF,0.2),quantile(all.d$BF, 0.8)))

all.80m <- all.d %>% group_by(variable) %>% 
  filter(between(GL,quantile(all.d$GL,0.1),quantile(all.d$GL, 0.9)) | 
           between(BF,quantile(all.d$BF,0.1),quantile(all.d$BF, 0.9)))


# Output Linear Regression Equation for all Hypotheses ------------------------------------------------------------------------------------------------------------------------------------------------------

# See Calcs file 

# Graph Function------------------------------------------------------------------------------------------------------------------------------------------------------

SL <- function(df,x.val,y.val,cc,lowerx,upperx,ticx,lowery,uppery,ticy){
  x.tic <- seq(-lowerx,upperx,ticx) 
  y.tic <- seq(-lowery,uppery,ticy)
  #x.tic <- c(-lowerx,upperx,ticx) 
  #y.tic <- c(-lowery,uppery,ticy)
  scat <- ggplot(df, aes(x=x.val,y=y.val)) + 
    geom_point(alpha=0.75, aes(color=as.factor(.data[['diff']])), size=0.75) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=12, color="black"),
      text = element_text(size=14),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5),
      legend.background = element_rect(colour = "transparent",fill = "transparent"),
      axis.title.x=element_blank(),
      axis.title.y=element_blank()) +
    #coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    #labs(x="dGLS",y="ln(BF)", color="Metric Disagreement") + 
    scale_color_manual(values=cc) 
  return(scat)
}

# Linear Regression Functions------------------------------------------------------------------------------------------------------------------------------------------------------

lm_eqn = function(xvals,yvals){
  m = lm(yvals ~ xvals) #y ~ x
  eq <- substitute(italic(y) == a + b %.% italic(x)*",  "~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2), 
                        b = format(unname(coef(m)[2]), digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

LR <- function(dd,h,cc,lowerx,upperx,ticx,lowery,uppery,ticy,xe,ye,ds){
  df <- dd %>% filter(variable == h)
  x.val <- df$GL
  y.val <- df$BF
  d.eqn <- data.frame(
    GL = xe,
    BF = ye,
    eqn = lm_eqn(df$GL,df$BF))
  
  lr <- SL(df,x.val,y.val,cc,lowerx,upperx,ticx,lowery,uppery,ticy) +
    geom_smooth(method=lm ,color="black", formula = y ~ x,size=0.25,se=TRUE) + 
    #geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.25) +
    #geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.25) + 
    theme(legend.position = "none") + ggtitle(ds) + 
    geom_text(data = d.eqn, aes(x = mean(range(x.val)), y = BF, label = eqn),color="black", 
              size=4, parse = TRUE)
  return(lr) 
}



# Linear Scatter TvS -----------------------------------------------------------------------------------------------------------------------------------------------------

  
tt <- paste(dataset,"Toxicofera v Scleroglossa",sep = " - ")
h <- 'TvS'

df <- all.d %>% filter(variable == h)
xval <- df$GL
yval <- df$BF
max(abs(xval))
max(abs(yval))

cc1 <- c("springgreen4","springgreen4","springgreen4","springgreen4","springgreen4")

quartz()
# x y
a.lr <- LR(all.d,h,cc1,60,60,10,100,100,20,0,90,"a. TvS - All")
a.lr

a.lr.80m <- LR(all.80m,h,cc1,60,60,5,100,100,20,0,80,"b. TvS - 80m") 
a.lr.80m

a.lr.50m <- LR(all.50m,h,cc1,60,60,5,100,100,10,0,45,"50m") 
#a.lr.50m

# named 2x in paper cuz its 2x threshold vut its +1 threshold above threshold.
a.lr.2x <- LR(all.1x,h,cc1,60,60,5,100,100,10,0,45,"2x") 
#a.lr.2x

a.lr.1x <- LR(all.0x,h,cc1,60,60,5,100,100,5,0,15,"1x") 
#a.lr.1x


a.lr.50l <- LR(all.50l,h,cc1,60,60,10,100,100,20,0,45,"50l") 
a.lr.50h <- LR(all.50h,h,cc1,60,60,10,100,100,20,0,90,"50h") 

f <- ggarrange(a.lr.50m,a.lr.50l,a.lr.50h, ncol=1, nrow=3,align="h")


lr <- ggarrange(a.lr,a.lr.80m,a.lr.2x,a.lr.1x,  ncol=2, nrow=2,align="h")


fig <- annotate_figure(lr,
                       bottom = text_grob("dGLS", size = 16),
                       left = text_grob("ln(BF)",  size = 16, rot = 90),
                       top = text_grob(tt,  size = 18))


fig50 <- annotate_figure(f,
                         bottom = text_grob("dGLS", size = 16),
                         left = text_grob("ln(BF)",  size = 16, rot = 90),
                         top = text_grob(tt,  size = 18))

# for manuscript 

m <- ggarrange(a.lr,a.lr.80m,a.lr.50m,a.lr.2x,  ncol=2, nrow=2,align="h")


figm <- annotate_figure(m,
                        bottom = text_grob("dGLS", size = 16),
                        left = text_grob("ln(BF)",  size = 16, rot = 90),
                        top = text_grob(tt,  size = 18))

m2 <- ggarrange(a.lr,a.lr.80m, ncol=2, nrow=1,align="h",labels = c("A", "B"))


m2


#ggsave(paste(dataset,h,"scatter_lrm.pdf",sep="_"), plot=figm,width = 8, height = 8, units = "in", device = 'pdf',bg = "transparent")


#ggsave(paste(dataset,h,"scatter_lr.pdf",sep="_"), plot=fig,width = 8, height = 8, units = "in", device = 'pdf',bg = "transparent")

#ggsave(paste(dataset,h,"scatter_lr50.pdf",sep="_"), plot=fig50,width = 4, height = 9, units = "in", device = 'pdf',bg = "transparent")



# Linear Scatter Tox -----------------------------------------------------------------------------------------------------------------------------------------------------



tt <- paste(dataset,"AI v SA",sep = " - ")
h <- 'AIvSA'
#tt <- paste(dataset,"AI v SI",sep = " - ")
#h <- 'AIvSI'
#tt <- paste(dataset,"SA v SI",sep = " - ")
#h <- 'SAvSI'

df <- all.d %>% filter(variable == h)
xval <- df$GL
yval <- df$BF
max(abs(xval))
max(abs(yval))

cc1 <- c("springgreen4","springgreen4","springgreen4","springgreen4","springgreen4")

quartz()
a.lr <- LR(all.d,h,cc1,40,40,2,60,60,5,0,15,"c. AIvSA - All")
a.lr
 
a.lr.80m <- LR(all.80m,h,cc1,40,40,2,60,60,5,0,15,"d. AIvSA - 80m") 
a.lr.80m

a.lr.50m <- LR(all.50m,h,cc1,40,40,2,60,60,5,0,10,"50m") 
#a.lr.50m

# named 2x in paper cuz its 2x threshold vut its +1 threshold above threshold.
a.lr.2x <- LR(all.1x,h,cc1,40,40,2,60,60,10,0,25,"2x") 
#a.lr.2x

a.lr.1x <- LR(all.0x,h,cc1,40,40,2,60,60,5,0,15,"1x") 
#a.lr.1x


a.lr.50l <- LR(all.50l,h,cc1,40,40,2,60,60,5,0,5,"50l") 
#a.lr.50l
a.lr.50h <- LR(all.50h,h,cc1,40,40,2,60,60,5,0,15,"50h") 

f <- ggarrange(a.lr.50m,a.lr.50l,a.lr.50h, ncol=1, nrow=3,align="h")
#f

lr <- ggarrange(a.lr,a.lr.80m,a.lr.2x,a.lr.1x,  ncol=2, nrow=2,align="h")
#lr



fig <- annotate_figure(lr,
                       bottom = text_grob("dGLS", size = 16),
                       left = text_grob("ln(BF)",  size = 16, rot = 90),
                       top = text_grob(tt,  size = 18))


fig50 <- annotate_figure(f,
                       bottom = text_grob("dGLS", size = 16),
                       left = text_grob("ln(BF)",  size = 16, rot = 90),
                       top = text_grob(tt,  size = 18))



# for manuscript 
m <- ggarrange(a.lr,a.lr.80m,a.lr.50m,a.lr.2x,  ncol=2, nrow=2,align="h")

figm <- annotate_figure(m,
                        bottom = text_grob("dGLS", size = 16),
                        left = text_grob("ln(BF)",  size = 16, rot = 90),
                        top = text_grob(tt,  size = 18))


m2t <- ggarrange(a.lr,a.lr.80m, ncol=2, nrow=1,align="h")


m2t2 <- ggarrange(m2,m2t,ncol=1, nrow=2)
m2t2

figm <- annotate_figure(m2t2,
                        bottom = text_grob("dGLS", size = 16),
                        left = text_grob("ln(BF)",  size = 16, rot = 90),
                        top = text_grob("Reeder",  size = 18))
figm


ggsave(paste(dataset,h,"scatter_lrm2.pdf",sep="_"), plot=figm,width = 8, height = 8, units = "in", device = 'pdf',bg = "transparent")


#ggsave(paste(dataset,h,"scatter_lr.pdf",sep="_"), plot=fig,width = 8, height = 8, units = "in", device = 'pdf',bg = "transparent")

#ggsave(paste(dataset,h,"scatter_lr50.pdf",sep="_"), plot=fig50,width = 4, height = 9, units = "in", device = 'pdf',bg = "transparent")



# Output disagreement summaries -------------------------------------------------Classify discordance----------------------------------------------------------------------------------------------------

a <- all %>% dplyr::rename(Hypothesis=variable) %>% 
  mutate(sig = case_when(between(BF,0,5) & between(GL,0,0.5) ~ 'neu.s',
                                    between(BF,-5,0) & between(GL,-0.5,0) ~ 'neu.s',
                                    between(BF,0,5) & between(GL,-0.5,0) ~ 'neu.o',
                                    between(BF,-5,0) & between(GL,0,0.5) ~ 'neu.o',
                                    BF >= 5 & GL >= 0.5 ~ 'sig.s',
                                    BF <= -5 & GL <= -0.5 ~ 'sig.s',
                                    BF <= -5 & GL >= 0.5 ~ 'sig.o',
                                    BF >= 5 & GL <= -0.5 ~ 'sig.o',
                                    BF <= -5 & between(GL,-0.5,0) ~ 'neuG.sigB.s',
                                    BF >= 5 & between(GL,0,0.5) ~ 'neuG.sigB.s',
                                    GL <= -0.5 & between(BF,-5,0) ~ 'neuB.sigG.s',
                                    GL >= 0.5 & between(BF,0,5) ~ 'neuB.sigG.s',
                                    BF <= -5 & between(GL,0,0.5) ~ 'neuG.sigB.o',
                                    BF >= 5 & between(GL,-0.5,0) ~ 'neuG.sigB.o',
                                    GL <= -0.5 & between(BF,0,5) ~ 'neuB.sigG.o',
                                    GL >= 0.5 & between(BF,-5,0) ~ 'neuB.sigG.o')) 
  
loci <- length((a %>% filter(Hypothesis == "TvS"))$Locus)
# Calculate metrics for all loci, all loci is divided by 4 times loci to account for all hypotheses
x <- a %>% group_by(sig) %>% dplyr::summarise(Hypothesis='ALL',n.loci=n()) %>%
  mutate(percent.loci =  round((n.loci/(loci*4))*100,digits=2))
y <- a %>% group_by(sig,Hypothesis) %>% dplyr::summarise(n.loci=n()) %>% 
  arrange(sig,Hypothesis) %>% 
  mutate(percent.loci= round((n.loci/loci)*100,digits = 2))
# This will produce an error because it removes factors, this is OKAY.
z <- bind_rows(x,y)


# Summarize how many loci have some form of conflict 
q.a <- a %>% 
  group_by(Locus,sig) %>% 
  dplyr::summarize(sig.diff.all=n_distinct(Hypothesis)) 

q.t <- a %>% filter(Hypothesis != "TvS") %>%
  group_by(Locus,sig) %>% 
  dplyr::summarize(sig.diff.tox=n_distinct(Hypothesis))

q.s <- a %>% filter(Hypothesis == "TvS") %>%
  group_by(Locus,sig) %>% 
  dplyr::summarize(sig.diff.tvs=n_distinct(Hypothesis)) 


q <- full_join(full_join(q.a,q.t,by="Locus"),q.s, by="Locus")

qq <- q %>% group_by(sig.diff.all,sig.diff.tox,sig.diff.tvs,sig) %>% 
  dplyr::summarize(n.loci=n()) %>%
  mutate(percent.loci =  round((n.loci/(loci*4))*100,digits=2)) %>% 
  filter(!sig %in% c('neu.s','sig.s'))

# output excel sheet

library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_significanceConflicts.xlsx",sep='')

write.xlsx2(as.data.frame(z), file=fname, sheetName="bySignificance", row.names=FALSE)
write.xlsx2(as.data.frame(qq), file=fname, sheetName="byHypothesis", append=TRUE, row.names=FALSE)


# neu.s both neutral in same direction
# neu.o both neutral in opposite directions
# sig.s both sig in same direction 
# sig.o both sig in opposite directions
# neuG.sigB.s dgls neutral, bf significant - both same direction
# neuB.sigG.s bf neutral, dgls significant - both same direction 
# neuG.sigB.o dgls neutral, bf significant - in opposite directions
# neuB.sigG.o bf neutral, dgls significant - in opposite directions
# can cluster by same or opposite direction if needed. 
# Disagree Scatters --------------------------------------------------------------------------------------------------------

#tt <- paste(dataset,"Toxicofera v Scleroglossa",sep = " - ")
#h <- 'TvS'
#tt <- paste(dataset,"AI v SA",sep = " - ")
#h <- 'AIvSA'
#tt <- paste(dataset,"AI v SI",sep = " - ")
#h <- 'AIvSI'
tt <- paste(dataset,"SA v SI",sep = " - ")
h <- 'SAvSI'

df <- all.d %>% filter(variable == h)
xval <- df$GL
yval <- df$BF

max(abs(xval))
max(abs(yval))


#quartz()
cc1 <- c('grey80','orange3','slateblue','black')
MD <- SL(df,xval,yval,cc1,15,15,5,25,25,10) + 
  geom_vline(xintercept=c(0.5,-0.5),color=c("black"), linetype="dashed", size=0.25) +
  geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.25) + 
  ggtitle(tt)

MD
MD.z <- SL(df,xval,yval,cc1,10,10,5,10,10,5) + 
  geom_vline(xintercept=c(0.5,-0.5),color=c("black"), linetype="dashed", size=0.25) +
  geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.25) + 
  ggtitle(tt)

MD.z

md <- ggarrange(MD,MD.z, ncol=2, nrow=1,align="v")

md

#ggsave(paste(dataset,h,"scatter_disagree.pdf",sep="_"), plot=md,width = 16, height = 6, units = "in", device = 'pdf',bg = "transparent")
# Reeder, already zoomed in
#ggsave(paste(dataset,h,"scatter_disagree.pdf",sep="_"), plot=MD,width = 8, height = 6, units = "in", device = 'pdf',bg = "transparent")




