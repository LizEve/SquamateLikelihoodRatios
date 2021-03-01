rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library("viridis")
# set dataset 

#dataset <- "SinghalOG"
dataset <- "Streicher"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Colors 
# https://www.w3schools.com/colors/colors_picker.asp
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS"))
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))

# Get total number of loci 
loci <- length(mLGL$Locus)

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10 - not including average calcs
#names(mLBF[7:10])
#mLBF[7:10] <- mLBF[7:10]/2
# for singhal OG 
names(mLBF[4])
mLBF[4] <- mLBF[4]/2


# Want to stack datasets so GL and BF columns have both comparisons for each hypothesis. 

reformatdf <- function(df,h,h1,s1,keepCols){
  # keepCols is a vector of the locus and support type column indices
  # Grab column number  
  c1 <- match(h1,names(df))
  # Grab and rename column, add hypothesis column
  a.df <- df[,append(keepCols,c1)] 
  names(a.df)[names(a.df) == h1] <- h
  a.df$Hypothesis <- rep(h1,length(a.df[,1]))
  c <- match(h,names(a.df))
  # Adjust direction of support if needed. ie if column is AIvSA, but you want to know SAvAI
  a.df[,c] <- a.df[,c]*(s1)
  return(a.df)
}

# reformat 
AI.g <- bind_rows(reformatdf(mLGL,"AI","AIvSA",1,c(1,13)),
                  reformatdf(mLGL,"AI","AIvSI",1,c(1,13)))

SA.g <- bind_rows(reformatdf(mLGL,"SA","AIvSA",-1,c(1,13)),
                  reformatdf(mLGL,"SA","SAvSI",1,c(1,13)))

SI.g <- bind_rows(reformatdf(mLGL,"SI","AIvSI",-1,c(1,13)),
                  reformatdf(mLGL,"SI","SAvSI",-1,c(1,13)))

AI.b <- bind_rows(reformatdf(mLBF,"AI","AIvSA",1,c(1,13)),
                  reformatdf(mLBF,"AI","AIvSI",1,c(1,13)))

SA.b <- bind_rows(reformatdf(mLBF,"SA","AIvSA",-1,c(1,13)),
                  reformatdf(mLBF,"SA","SAvSI",1,c(1,13)))

SI.b <- bind_rows(reformatdf(mLBF,"SI","AIvSI",-1,c(1,13)),
                  reformatdf(mLBF,"SI","SAvSI",-1,c(1,13)))

TS.g <- reformatdf(mLGL,"TS","TvS",1,c(1,13)) %>% select(-SIvAvg)
TS.b <- reformatdf(mLBF,"TS","TvS",1,c(1,13)) %>% select(-SIvAvg)

# Singhal OG
#TS.g <- reformatdf(mLGL,"TS","TvS",1,c(1,3))
#TS.b <- reformatdf(mLBF,"TS","TvS",1,c(1,3))

# MERGE 

AI.a <- left_join(AI.g, AI.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
SA.a <- left_join(SA.g, SA.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
SI.a <- left_join(SI.g, SI.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))
TS.a <- left_join(TS.g, TS.b, by=c("Locus","Hypothesis"), suffix = c(".g", ".b"))

rm(AI.g,AI.b,SA.b,SA.g,SI.b,SI.g,TS.b,TS.g)

library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_TvSLoci.xlsx",sep='')

#write.xlsx2(as.data.frame(TS.a), file=fname, sheetName="loci", row.names=FALSE)



#------------------------------------------------------------------------------------------------------------------------------------------------------

# Quantiles ---------------------------------------------- Scatter Quantiles  --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Calcs------------------------------------------------------------------------------------------------------------------------------------------------------

# Adding two rows that give the difference in percentiles and rank between bf and gl 
# This will give the same rank number if the value is identicle in two cells - 
# hint - decile scatter is pointless, deleted from script
# Rank difference needs to be calulated differently because the total rank spots for df and gl are different
# Normalizing rank by calculating max rank by hypothesis - rank.g/max(rank.g) - using these fractions to compare bf and gl 
# Now rank difference is the percent difference? lets just call it normalized difference


AI.p <- AI.a %>%
  group_by(Hypothesis) %>% 
  mutate(per.g = ntile(desc(AI.g),100)) %>% 
  mutate(rank.g = dense_rank(desc(AI.g))) %>%
  mutate(per.b = ntile(desc(AI.b),100)) %>% 
  mutate(rank.b =dense_rank(desc(AI.b))) %>%
  mutate(per.diff = abs(per.g-per.b)) %>%
  mutate(rank.diff = abs((rank.g/max(rank.g))-(rank.b/max(rank.b)))*100)



SA.p <- SA.a %>%
  group_by(Hypothesis) %>% 
  mutate(per.g = ntile(desc(SA.g),100)) %>% 
  mutate(rank.g =dense_rank(desc(SA.g))) %>%
  mutate(per.b = ntile(desc(SA.b),100)) %>% 
  mutate(rank.b =dense_rank(desc(SA.b)))   %>%
  mutate(per.diff = abs(per.g-per.b)) %>%
  mutate(rank.diff = abs((rank.g/max(rank.g))-(rank.b/max(rank.b)))*100)



SI.p <- SI.a %>%
  group_by(Hypothesis) %>% 
  mutate(per.g = ntile(desc(SI.g),100)) %>% 
  mutate(rank.g = dense_rank(desc(SI.g))) %>%
  mutate(per.b = ntile(desc(SI.b),100)) %>% 
  mutate(rank.b = dense_rank(desc(SI.b)))   %>%
  mutate(per.diff = abs(per.g-per.b))  %>%
  mutate(rank.diff = abs((rank.g/max(rank.g))-(rank.b/max(rank.b)))*100)



TS.p <- TS.a %>%
  group_by(Hypothesis) %>% 
  mutate(per.g = ntile(desc(TS.g),100)) %>% 
  mutate(rank.g =dense_rank(desc(TS.g))) %>%
  mutate(per.b = ntile(desc(TS.b),100)) %>% 
  mutate(rank.b =dense_rank(desc(TS.b)))  %>%
  mutate(per.diff = abs(per.g-per.b)) %>%
  mutate(rank.diff = abs((rank.g/max(rank.g))-(rank.b/max(rank.b)))*100)


# Graph Function ------------------------------------------------------------------------------------------------------------------------------------------------------

SC <- function(df,x.tic,y.tic,legtit,l){
  # defining columns here to avoid repeating this over and over
  xcol <- match(paste(l,".g",sep = ''),names(df))
  ycol <- match(paste(l,".b",sep = ''),names(df))
  CC <- match(paste(l,".diff",sep = ''),names(df))
  
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]])) + 
    geom_point(alpha=1, aes(color=df[[CC]]), size=1) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=6, color="black"),
      text = element_text(size=10),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x=paste('dGLS ',l, sep=''),y=paste('ln(BF) ',l,sep=''),color=legtit) + 
    guides(colour = guide_legend(override.aes = list(size=2))) +
    scale_color_viridis(limits=c(0,100),direction = -1,option='D')
  #scale_color_gradientn(colors=heat.colors(2),limits=c(0,100))
  #scale_color_gradient2(low="blue", high="red",midpoint=mid)
  
  return(scat)
}


# All Rank #------------------------------------------------------------------------------------------------------------------------------------------------------

# BF vs dGLS colored by normalized difference in rank between dGLS and BF
# Each comparison has a graph. AIvSA, AIvSI, SAvSI, TvS

# Get scale for x and y axes
df1 <- AI.p[AI.p$Hypothesis=='AIvSA',] 
df2 <- AI.p[AI.p$Hypothesis=='AIvSI',] 
df3 <- SA.p[SA.p$Hypothesis=='SAvSI',] 
df4 <- TS.p

max(c(df2$rank.b,df2$rank.g))
max(c(df3$rank.b,df3$rank.g))
max(c(df4$rank.b,df4$rank.g))


# Tic marks 
t.rx <- seq(0,3000,500)
t.ry <- seq(0,3000,500)
ts.rx <- seq(0,5000,500)
ts.ry <- seq(0,5000,500)

#values,tick marks,column for values to color, legend title, rank/decile for axis labels
p1.a <- SC(df1,t.rx,t.ry,"normalized \n difference \n in rank","rank") + ggtitle(df1$Hypothesis)
p2.a <- SC(df2,t.rx,t.ry,"normalized \n difference \n  in rank","rank") + ggtitle(df2$Hypothesis)
p3.a <- SC(df3,t.rx,t.ry,"normalized \n difference \n in rank","rank") + ggtitle(df3$Hypothesis)
p4.a <- SC(df4,ts.rx,ts.ry,"normalized \n difference \n in rank","rank") + ggtitle(df4$Hypothesis)

p.a <- ggarrange(p1.a, p2.a, p3.a, p4.a, ncol=2, nrow=2, align="h")
#quartz()
p.a

#ggsave(paste(dataset,"_scatter_rank_diff.pdf",sep=""), plot=p.a,width = 9, height = 6, units = "in", device = 'pdf',bg = "transparent")
#ggsave(paste(dataset,"_scatter_rank_diff.pdf",sep=""), plot=p4.a,width = 9, height = 6, units = "in", device = 'pdf',bg = "transparent")

# All Percentiles #------------------------------------------------------------------------------------------------------------------------------------------------------

# BF vs dGLS colored by normalized difference in rank between dGLS and BF
# Each comparison has a graph. AIvSA, AIvSI, SAvSI, TvS


# Tic marks 
t.rx <- seq(0,100,20)
t.ry <- seq(0,100,20)
ts.rx <- seq(0,100,20)
ts.ry <- seq(0,100,20)

#values,tick marks,column for values to color, legend title, per/decile for axis labels
r1.a <- SC(df1,t.rx,t.ry,"difference in \n percentile","per") + ggtitle(df1$Hypothesis)
r2.a <- SC(df2,t.rx,t.ry,"difference in \n percentile","per") + ggtitle(df2$Hypothesis)
r3.a <- SC(df3,t.rx,t.ry,"difference in \n percentile","per") + ggtitle(df3$Hypothesis)
r4.a <- SC(df4,ts.rx,ts.ry,"difference in \n percentile","per") + ggtitle(df4$Hypothesis)

r.a <- ggarrange(r1.a, r2.a, r3.a, r4.a, ncol=2, nrow=2, align="h")
#quartz()
r.a

#ggsave(paste(dataset,"_scatter_percentile_diff.pdf",sep=""), plot=r.a,width = 9, height = 6, units = "in", device = 'pdf',bg = "transparent")
#ggsave(paste(dataset,"_scatter_percentile_diff.pdf",sep=""), plot=r4.a,width = 9, height = 6, units = "in", device = 'pdf',bg = "transparent")

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Compare hypotheses ------------------------------------ Loci across hypotheses --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Compare hypotheses to see if same loci give same support across hypotheses. ie, if it supports AI in AIvSI will it support it in AIvSA

# Graph Function------------------------------------------------------------------------------------------------------------------------------------------------------

SL <- function(df,x.val,y.val,x.tic,y.tic,cc,x.y){
  scat <- ggplot(df, aes(x=x.val,y=y.val)) + 
    geom_point(alpha=1, aes(color=as.factor(df[[x.y[3]]])), size=0.5) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=10, color="black"),
      text = element_text(size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5),
      legend.position = 'none'
    ) +
    coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x=x.y[1],y=x.y[2]) + 
    scale_color_manual(values=cc) 
  return(scat)
}

# dGLS-------------------------------------------------------- Loci across hypotheses ----------------------------------------------------------------------------------------------

# Set up hypotheses with correct directionality and add codes for one/both/neither with support for hypoth


GL <- mLGL %>% 
  select(Locus,TvS,AIvSA,AIvSI,SAvSI) %>% 
  mutate(SAvAI=AIvSA*-1) %>%
  mutate(SIvAI=AIvSI*-1) %>%
  mutate(SIvSA=SAvSI*-1) %>%
  mutate(AI = case_when(AIvSI < 0.5 &  AIvSA < 0.5 ~ 0,
                        AIvSI >= 0.5 & AIvSA >= 0.5 ~ 2,
                        AIvSI >= 0.5 & AIvSA < 0.5 ~ 1,
                        AIvSI < 0.5 & AIvSA >= 0.5 ~ 1)) %>%
  mutate(SA = case_when(SAvSI < 0.5 &  SAvAI < 0.5 ~ 0,
                        SAvSI >= 0.5 & SAvAI >= 0.5 ~ 2,
                        SAvSI >= 0.5 & SAvAI < 0.5 ~ 1,
                        SAvSI < 0.5 & SAvAI >= 0.5 ~ 1)) %>%
  mutate(SI = case_when(SIvSA < 0.5 &  SIvAI < 0.5 ~ 0,
                        SIvSA >= 0.5 & SIvAI >= 0.5 ~ 2,
                        SIvSA >= 0.5 & SIvAI < 0.5 ~ 1,
                        SIvSA < 0.5 & SIvAI >= 0.5 ~ 1))

# Graph 

df <- GL
max(c(abs(df$AIvSA),abs(df$AIvSI),abs(df$SAvSI)))
lower <- 170 # 45
upper <- 170
x.t <- seq(-lower,upper,50) # 15
y.t <- seq(-lower,upper,50)

df <- GL
x.v1 <- GL$AIvSA
y.v1 <- GL$AIvSI
x.y1 <- c("AIvSA","AIvSI","AI")
cc1 <- c('grey56','#1e7b59','#46d29f')
quartz()
ai.g <- SL(df,x.v1,y.v1,x.t,y.t,cc1,x.y1) + geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.25) +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.25)

x.v2 <- GL$SAvAI
y.v2 <- GL$SAvSI
x.y2 <- c("SAvAI","SAvSI","SA")
cc2 <- c('grey56','#243a5b','#4975b6')
sa.g <- SL(df,x.v2,y.v2,x.t,y.t,cc2,x.y2) + geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.25) +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.25)


x.v3 <- GL$SIvAI
y.v3 <- GL$SIvSA
x.y3 <- c("SIvAI","SIvSA","SI")
cc3 <- c('grey56','#4d4d00','#cccc00')
si.g <- SL(df,x.v3,y.v3,x.t,y.t,cc3,x.y3) + geom_hline(yintercept=c(0.5),color=c("black"), linetype="dashed", size=0.25) +
  geom_vline(xintercept=c(0.5),color=c("black"), linetype="dashed", size=0.25)

g <- annotate_figure(ggarrange(ai.g,sa.g,si.g, ncol=1, nrow=3, align="v",labels="auto"),
                     top = text_grob("dGLS", color = "black", size = 16))
g

# BF --------------------------------------------------------- Loci across hypotheses ---------------------------------------------------------------------------------------------

# Set up hypotheses with correct directionality and add codes for one/both/neither with support for hypoth

BF <- mLBF %>% 
  select(Locus,TvS,AIvSA,AIvSI,SAvSI) %>% 
  mutate(SAvAI=AIvSA*-1) %>%
  mutate(SIvAI=AIvSI*-1) %>%
  mutate(SIvSA=SAvSI*-1) %>%
  mutate(AI = case_when(AIvSI < 5 &  AIvSA < 5 ~ 0,
                        AIvSI >= 5 & AIvSA >= 5 ~ 2,
                        AIvSI >= 5 & AIvSA < 5 ~ 1,
                        AIvSI < 5 & AIvSA >= 5 ~ 1)) %>%
  mutate(SA = case_when(SAvSI < 5 &  SAvAI < 5 ~ 0,
                        SAvSI >= 5 & SAvAI >= 5 ~ 2,
                        SAvSI >= 5 & SAvAI < 5 ~ 1,
                        SAvSI < 5 & SAvAI >= 5 ~ 1)) %>%
  mutate(SI = case_when(SIvSA < 5 &  SIvAI < 5 ~ 0,
                        SIvSA >= 5 & SIvAI >= 5 ~ 2,
                        SIvSA >= 5 & SIvAI < 5 ~ 1,
                        SIvSA < 5 & SIvAI >= 5 ~ 1))


# Graph

max(c(abs(BF$AIvSA),abs(BF$AIvSI),abs(BF$SAvSI)))
lower <- 90 
upper <- 90
x.t <- seq(-lower,upper,20) 
y.t <- seq(-lower,upper,20)

df2 <- BF
x.v4 <- BF$AIvSA
y.v4 <- BF$AIvSI
x.y4 <- c("AIvSA","AIvSI","AI")
cc4 <- c('grey56','#1e7b59','#46d29f')
#quartz()
ai.b <- SL(df2,x.v4,y.v4,x.t,y.t,cc4,x.y4) + geom_hline(yintercept=c(5),color=c("black"), linetype="dashed", size=0.25) +
  geom_vline(xintercept=c(5),color=c("black"), linetype="dashed", size=0.25)
ai.b

x.v5 <- BF$SAvAI
y.v5 <- BF$SAvSI
x.y5 <- c("SAvAI","SAvSI","SA")
cc5 <- c('grey56','#243a5b','#4975b6')
sa.b <- SL(df2,x.v5,y.v5,x.t,y.t,cc5,x.y5) + geom_hline(yintercept=c(5),color=c("black"), linetype="dashed", size=0.25) +
  geom_vline(xintercept=c(5),color=c("black"), linetype="dashed", size=0.25)


x.v6 <- BF$SIvAI
y.v6 <- BF$SIvSA
x.y6 <- c("SIvAI","SIvSA","SI")
cc6 <- c('grey56','#4d4d00','#cccc00')
si.b <- SL(df2,x.v6,y.v6,x.t,y.t,cc6,x.y6) + geom_hline(yintercept=c(5),color=c("black"), linetype="dashed", size=0.25) +
  geom_vline(xintercept=c(5),color=c("black"), linetype="dashed", size=0.25)

b <- annotate_figure(ggarrange(ai.b,sa.b,si.b, ncol=1, nrow=3, align="v",labels=c("d","e","f")),
                     top = text_grob("ln(BF)", color = "black", size = 16))
b


f <- ggarrange(g,b, ncol=2, nrow=1, align="v")
f

#ggsave(paste(dataset,"_scatter_loci_compare.pdf",sep=""), plot=f,width = 9, height = 12, units = "in", device = 'pdf',bg = "transparent")
#ggsave(paste(dataset,"_scatter_loci_comparezoom.pdf",sep=""), plot=f,width = 9, height = 12, units = "in", device = 'pdf',bg = "transparent")


#------------------------------------------------------------------------------------------------------------------------------------------------------

# CSV files ------------------------------------ CSV files --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------


# Smash data ----------------------------------------- Significant disagree CSV------------------------------------------------------------------------------------------------------------

# Look at "significance" for dGLS and BF. Record which loci conflict

# Create dataframe with all metadata, locus names, one col for hypothesis, one call for bf, and one gl
keep <- c(1,seq(14,22))

B <- bind_rows(reformatdf(mLBF,"BF","AIvSA",1,keep),
               reformatdf(mLBF,"BF","AIvSI",1,keep),
               reformatdf(mLBF,"BF","SAvSI",1,keep),
               reformatdf(mLBF,"BF","TvS",1,keep))

G <- bind_rows(reformatdf(mLGL,"GL","AIvSA",1,keep),
               reformatdf(mLGL,"GL","AIvSI",1,keep),
               reformatdf(mLGL,"GL","SAvSI",1,keep),
               reformatdf(mLGL,"GL","TvS",1,keep))
# for SinghalOG 
#keep <- c(1,seq(5,13)) 
#B <- reformatdf(mLBF,"BF","TvS",1,keep)
#G <- reformatdf(mLGL,"GL","TvS",1,keep)

metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")
all <- left_join(B, G, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))


# Classify discordance ----------------------------------------- Significant disagree CSV------------------------------------------------------------------------------------------------------------

# neu.s both neutral in same direction
# neu.o both neutral in opposite directions
# sig.s both sig in same direction 
# sig.o both sig in opposite directions
# neuG.sigB.s dgls neutral, bf significant - both same direction
# neuB.sigG.s bf neutral, dgls significant - both same direction 
# neuG.sigB.o dgls neutral, bf significant - in opposite directions
# neuB.sigG.o bf neutral, dgls significant - in opposite directions
# can cluster by same or opposite direction if needed. 


a <- all %>% mutate(sig = case_when(between(BF,0,5) & between(GL,0,0.5) ~ 'neu.s',
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


# Classify discordance ----------------------------------------- Rank/Percentile disagree CSV------------------------------------------------------------------------------------------------------------

# calculate percentage for whole dataset
# add normalization to rank difference calculation

p <- a %>%
  group_by(Hypothesis) %>% 
  mutate(per.g = ntile(desc(GL),100)) %>% 
  mutate(rank.g = dense_rank(desc(GL))) %>%
  mutate(per.b = ntile(desc(BF),100)) %>% 
  mutate(rank.b = dense_rank(desc(BF))) %>%
  mutate(per.diff = (per.g-per.b)) %>%
  mutate(rank.diff = ((rank.g/max(rank.g))-(rank.b/max(rank.b)))*100)

# Rank diff - g - b. positive - higher G. negative - higher B
# same - within 50% of rank/percentile 
# rank.b - b is more than 50% greater than g in terms of rank 
# rank.g - g is more than 50% greater than b in terms of rank
# per.b - b is more than 50 percentiles greater than g 
# per.g - g is more than 50 percentiles greater than b 


# Add columns for whether normalized b and g rank or percentile difference is larger than half the total ranks
p <- p %>% mutate(rank = case_when(rank.diff >= 50 ~ 'rank.g',
                                   rank.diff <= -50 ~ 'rank.b',
                                   between(rank.diff,-50,50) ~ 'same')) %>%
  mutate(per = case_when(per.diff >= 50 ~ 'per.g',
                         per.diff <= -50 ~ 'per.b',
                         between(per.diff,-50,50) ~ 'same'))

# Summarize number of loci for each type of discordance
s.sum <- bind_rows(
  p %>% group_by(sig) %>% dplyr::summarise(Hypothesis='ALL',n.loci=n()),
  p %>% group_by(sig,Hypothesis) %>% dplyr::summarise(n.loci=n()) %>% 
    arrange(sig,Hypothesis)) %>% 
  mutate(prop.loci= case_when(Hypothesis == "ALL" ~ round((n.loci/(loci*4))*100,digits=2),
                              Hypothesis != "ALL" ~ round((n.loci/loci)*100,digits = 2)))

r.sum <- bind_rows(
  p %>% group_by(rank) %>% summarise(Hypothesis='ALL',n.loci=n()),
  p %>% group_by(rank,Hypothesis) %>% summarise(n.loci=n()) %>% 
    arrange(rank,Hypothesis)) %>% 
  mutate(prop.loci= case_when(Hypothesis == "ALL" ~ round((n.loci/(loci*4))*100,digits=2),
                              Hypothesis != "ALL" ~ round((n.loci/loci)*100,digits = 2)))

p.sum <- bind_rows(
  p %>% group_by(per) %>% summarise(Hypothesis='ALL',n.loci=n()),
  p %>% group_by(per,Hypothesis) %>% summarise(n.loci=n()) %>% 
    arrange(per,Hypothesis)) %>% 
  mutate(prop.loci= case_when(Hypothesis == "ALL" ~ round((n.loci/(loci*4))*100,digits=2),
                              Hypothesis != "ALL" ~ round((n.loci/loci)*100,digits = 2)))


# Summarize ----------------------------------------- Disagree CSV------------------------------------------------------------------------------------------------------------


# List loci and the number of hypotheses they have issues across. ie issue for all hypotheses or only one. 
y1 <- p %>% 
  filter(!sig %in% c('neu.s','sig.s')) %>% 
  group_by(Locus) %>% 
  summarize(n.h.sig=n_distinct(Hypothesis))

y2 <- p %>% 
  filter(!rank %in% c('same')) %>% 
  group_by(Locus) %>% 
  summarize(n.h.rank=n_distinct(Hypothesis))


y3 <- p %>% 
  filter(!per %in% c('same')) %>% 
  group_by(Locus) %>% 
  summarize(n.h.per=n_distinct(Hypothesis))

y <- full_join(full_join(y1,y2,by='Locus'),y3,by='Locus')


# Output ----------------------------------------- Disagree CSV------------------------------------------------------------------------------------------------------------

# Want to out put list of loci - y, and summary tables - r,p,s.sum 

library(xlsx)

fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_disagreeingLociSummary.xlsx",sep='')

x <- full_join(p,y,by="Locus") %>% filter(rank.diff<=-10| rank.diff>=10)
x$rank.diff.abs <- abs(x$rank.diff)
x$Locus.2 <- x$Locus
options(java.parameters = "-Xmx1024m")
library(XLConnect)
write.xlsx2(as.data.frame(x), file=fname, sheetName="all", row.names=FALSE)
write.xlsx2(as.data.frame(y), file=fname, sheetName="loci", append=TRUE,row.names=FALSE)
write.xlsx2(as.data.frame(s.sum), file=fname, sheetName="sig", append=TRUE, row.names=FALSE)
write.xlsx2(as.data.frame(r.sum), file=fname, sheetName="rank", append=TRUE, row.names=FALSE)
write.xlsx2(as.data.frame(p.sum), file=fname, sheetName="per", append=TRUE, row.names=FALSE)


# old 
#write.csv(z,paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_disagreeingLoci.csv",sep=''))
#write.csv(y,paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_disLociSummary.csv",sep=''))



#------------------------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------- Disagree Scatters --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Graph Function ------------------------------------------------------------------------------------------------------------------------------------------------------

SC <- function(df,xcol,ycol,c,x.tic,y.tic,xl,yl){
  scat <- ggplot(df, aes(x=xcol,y=ycol)) + 
    geom_point(alpha=1, size=2, aes(color=c)) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=6, color="black"),
      text = element_text(size=10),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x=xl,y=yl,color="abs val of x") + 
    scale_color_viridis(limits=c(min(c),max(c)),direction = 1,option='D')
  #guides(colour = guide_legend(override.aes = list(size=2))) +
  #scale_color_viridis(limits=c(-100,100),direction = -1,option='D')
  #scale_color_gradientn(colors=heat.colors(2),limits=c(0,100))
  #scale_color_gradient2(low="blue", high="red",midpoint=mid)
  
  return(scat)
}


# Graph BF GL Scatters------------------------------------------------------------------------------------------------------------------------------------------------------

all <- p

# Rank
x <- abs(all$rank.diff)
c <- all$rank.diff
max(x)
min(x)
ts.rx <- seq(0,100,20)
xlab <- "normalized rank difference"

# Bf and GL #####

xlab <- "ln(BF)"
x <- all$BF
max(x)
min(x)
ts.rx <- seq(-100,60,20)


ylab <- "GL"
y <- all$GL
max(y)
min(y)
ts.ry <- seq(-80,30,10)
c <- abs(all$rank.diff)
r <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) + labs(color="abs rank diff")
c <- all$rank.diff
r2 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) + labs(color="rank diff")


z <- ggarrange(r,r2, ncol=1, nrow=2,align="v")
z

#ggsave(paste(dataset,"_rank_bfgl.pdf",sep=""), plot=z ,width = 4, height = 7, units = "in", device = 'pdf',bg = "transparent")

# Graph diff Histogram ---------------------------------------------- Histogram --------------------------------------------------------------------------------------------------------

H <- function(df,xval,cc,hbin,xlab){
  ggplot(data=df, aes(x=xval)) +
    geom_histogram(bins=hbin, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1) +
    #geom_histogram(breaks=brx, alpha=1, position="dodge",  fill=cc, color="grey",size=0.1)+ 
    #geom_histogram(binwidth = max(abs(xval))*0.01, alpha=1, position="dodge",color="grey", fill=cc, size=0.1)+ 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=10),
          axis.text = element_text(size=6, color="black"),
          text = element_text(size=8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)) +
    labs(y="Number of Loci",x=xlab)  + 
    #coord_cartesian(xlim = x.tic) +
    #coord_cartesian(ylim=ytic, xlim = xtic) +
    #scale_y_continuous() +
    #scale_x_continuous(breaks = xtic) +
    # geom_vline(xintercept=lines.g,color=c("black"), linetype="dashed", size=0.5) +
    geom_vline(xintercept=c(0),color=c("black"), linetype="dashed", size=0.2) 
}


#quartz()

h.bin <- 30
hy <- "TvS"
#hy <- "All"
#hy <- "AIvSA"
#all <- p[p$Hypothesis == "AIvSA",]
all <- p[p$Hypothesis == "TvS",]
#all <- p

# Set input data 

x.val <- all$rank.diff
df <- all
cc <- "#2BB07FFF"
# Get max min for graph to set x axis values
max(abs(x.val))
m <- 100
x.tic <- seq(-m,m,20)

H(df,x.val,cc,h.bin,'dGLS values') +
  coord_cartesian(xlim = x.tic) + 
  scale_x_continuous(breaks = x.tic) 

y.tic <- seq(0,900,100)

rh <- H(df,x.val,cc,h.bin,'Difference in rank. - BF higher + GL higher') +
  coord_cartesian(xlim = x.tic, ylim=y.tic) + 
  scale_x_continuous(breaks = x.tic) + 
  scale_y_continuous(breaks = y.tic)
rh

x.val <- all$per.diff
df <- all
cc <- "#2BB07FFF"
# Get max min for graph to set x axis values
#max(abs(x.val))
#m <- 70
#x.tic <- seq(-m,m,10)

H(df,x.val,cc,h.bin,'dGLS values') +
  coord_cartesian(xlim = x.tic) + 
  scale_x_continuous(breaks = x.tic) 

#y.tic <- seq(0,40,5)

ph <- H(df,x.val,cc,h.bin,'Difference in percentile. - BF higher + GL higher') +
  coord_cartesian(xlim = x.tic, ylim=y.tic) + 
  scale_x_continuous(breaks = x.tic) + 
  scale_y_continuous(breaks = y.tic)
ph


title <- paste(dataset," - ",hy,sep="") 

h <- annotate_figure(ggarrange(rh,ph, ncol=1, nrow=2,align="v"),
                     top = text_grob(title, color = "black", size = 16))

h

#ggsave(paste(dataset,"_rank_hist_",hy,".pdf",sep=""), plot=h ,width = 5, height = 4, units = "in", device = 'pdf',bg = "transparent")



# Graph Metadata Scatters------------------------------------------------------------------------------------------------------------------------------------------------------

#all <- p[p$Hypothesis == "AIvSA",]
all <- p[p$Hypothesis == "TvS",]
#all <- p


ylab <- "num taxa"
y <- all$GL
max(y)
min(y)
ts.ry <- seq(200,290,30)
gl <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
gl 

ylab <- "alignment score"
y <- all$MEAN_COL_SCORE
ts.ry <- seq(0,1,0.1)
g1 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
g1 

ylab <- "num taxa"
y <- all$Sequences
max(y)
min(y)
ts.ry <- seq(200,290,30)
g2 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
g2 

ylab <- "alignment length"
y <- all$Columns
max(y)
min(y)
ts.ry <- seq(300,2100,300)
g3 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
g3 

names(all)
ylab <- "PI"
y <- all$Pars_Info
max(y)
min(y)
ts.ry <- seq(100,1500,100)
g4 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
g4

names(all)
ylab <- "Singleton Sites"
y <- all$Sing_Sites
max(y)
min(y)
ts.ry <- seq(0,200,50)
g5 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
g5

names(all)
ylab <- "Cons_Sites"
y <- all$Cons_Sites
max(y)
min(y)
ts.ry <- seq(0,1200,100)
g6 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
g6

ylab <- "Chi2_Fail"
y <- all$Chi2_Fail
max(y)
min(y)
ts.ry <- seq(0,250,50)
g7 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
g7

ylab <- "Gaps"
y <- all$Gaps_Ambig
max(y)
min(y)
ts.ry <- seq(0,6,1)
g8 <- SC(all,x,y,c,ts.rx,ts.ry,xlab,ylab) 
g8

g <- ggarrange(g1,g2,g3,g4,g5,g6,g7,g8, ncol=2, nrow=4,align="v")
g


ggsave(paste(dataset,"_rank_metaData.pdf",sep=""), plot=g ,width = 6, height = 10, units = "in", device = 'pdf',bg = "transparent")



#------------------------------------------------------------------------------------------------------------------------------------------------------

# Testin ------------------------------------ Testing --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------------------------------------------

# Mean Mode Variance ------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

bf <- all %>% group_by(Hypothesis) %>%
  summarise(Mean=mean(BF), Median=median(BF), Std=sd(BF), Mad=mad(BF), IQR=IQR(BF), Max=max(BF), Min=min(BF))

gl <- all %>% group_by(Hypothesis) %>%
  summarise(Mean=mean(GL), Median=median(GL), Std=sd(GL), Mad=mad(GL), IQR=IQR(GL), Max=max(GL), Min=min(GL))

# copy and paste into metadata excel sheet for now. 

#Outliers based on IQR------------------------------------------------------------------------------------------------------------------------------------------------------

# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r
# https://www.r-bloggers.com/combined-outlier-detection-with-dplyr-and-ruler/


data <- all$BF
lowerq = quantile(data)[2]
upperq = quantile(data)[4]
iqr = upperq - lowerq 
extreme.threshold.upper = (iqr * 3) + upperq
extreme.threshold.lower = lowerq - (iqr * 3)

x <- all %>% group_by(Hypothesis) %>% mutate(bo = case_when(BF <= extreme.threshold.lower ~ TRUE,
                                                            BF >= extreme.threshold.upper ~ TRUE,
                                                            between(BF,extreme.threshold.lower,extreme.threshold.upper) ~ FALSE)) %>%
  count(bo)

y <- all %>% group_by(Hypothesis) %>% mutate(go = case_when(GL <= extreme.threshold.lower ~ TRUE,
                                                            GL >= extreme.threshold.upper ~ TRUE,
                                                            between(GL,extreme.threshold.lower,extreme.threshold.upper) ~ FALSE)) %>%
  count(go)



