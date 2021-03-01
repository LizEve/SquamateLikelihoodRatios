rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# set dataset 

dataset <- "Singhal"


setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")

# Add column to name type of support 
mLGL <- mutate(mLGL, supportType = case_when(TvS != 'a' ~ "dGLS"))
mLBF <- mutate(mLBF, supportType = case_when(TvS != 'a' ~ "BF"))

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10
mLBF[7:10] <- mLBF[7:10]/2

# Check number of NA 
sum(is.na(mLBF$MEAN_COL_SCORE))

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
keep <- c(1,seq(14,22))
AI.g <- bind_rows(reformatdf(mLGL,"AI","AIvSA",1,keep),
                  reformatdf(mLGL,"AI","AIvSI",1,keep))

SA.g <- bind_rows(reformatdf(mLGL,"SA","AIvSA",-1,keep),
                  reformatdf(mLGL,"SA","SAvSI",1,keep))

SI.g <- bind_rows(reformatdf(mLGL,"SI","AIvSI",-1,keep),
                  reformatdf(mLGL,"SI","SAvSI",-1,keep))

AI.b <- bind_rows(reformatdf(mLBF,"AI","AIvSA",1,keep),
                  reformatdf(mLBF,"AI","AIvSI",1,keep))

SA.b <- bind_rows(reformatdf(mLBF,"SA","AIvSA",-1,keep),
                  reformatdf(mLBF,"SA","SAvSI",1,keep))

SI.b <- bind_rows(reformatdf(mLBF,"SI","AIvSI",-1,keep),
                  reformatdf(mLBF,"SI","SAvSI",-1,keep))

TS.g <- reformatdf(mLGL,"TS","TvS",-1,keep)
TS.b <- reformatdf(mLBF,"TS","TvS",-1,keep)

# MERGE - yea i could have merged first, but the function was already written and i didnt want to mess with it. 

AI.a <- left_join(AI.g, AI.b, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))
SA.a <- left_join(SA.g, SA.b, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))
SI.a <- left_join(SI.g, SI.b, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))
TS.a <- left_join(TS.g, TS.b, by=c("Locus","Hypothesis",metaData), suffix = c(".g", ".b"))

sapply(AI.a,typeof)

# Smash all together, then smash bf and gl together each in 1 column
all <- bind_rows(bind_rows(bind_rows(AI.a,SA.a),SI.a),TS.a) %>% 
  mutate(BF=coalesce(AI.b,SA.b,SI.b,TS.b)) %>% 
  mutate(GL=coalesce(AI.g,SA.g,SI.g,TS.g)) %>%
  select(-c(AI.b,SA.b,SI.b,TS.b,AI.g,SA.g,SI.g,TS.g))

TX <- bind_rows(bind_rows(AI.a,SA.a),SI.a) %>% 
  mutate(BF=coalesce(AI.b,SA.b,SI.b)) %>% 
  mutate(GL=coalesce(AI.g,SA.g,SI.g)) %>%
  select(-c(AI.b,SA.b,SI.b,AI.g,SA.g,SI.g))


TS <- TS.a %>% 
  mutate(BF=coalesce(TS.b)) %>% 
  mutate(GL=coalesce(TS.g)) %>%
  select(-c(TS.b,TS.g))


# Check number of NA 
sum(is.na(all$BF))
sum(is.na(all$GL))

rm(AI.g,AI.b,SA.b,SA.g,SI.b,SI.g,TS.b,TS.g)

v <- 'tx'
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/",dataset,"_Corr_",v,"_rsvals.RData",sep=''))
#------------------------------------------------------------------------------------------------------------------------------------------------------

# Scatter ---------------------------------------------- Scatter BF vs dGLS --------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------

## Scatter graph, colored by hypothesis 

S <- function(df,xcol,ycol,x.tic,y.tic,cc,xlab,ylab){
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]])) + 
    geom_point(alpha=1, aes(color=df$Hypothesis), size=0.5) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=14, color="black"),
      text = element_text(size=14),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x=xlab,y=ylab,color="Hypothesis") + 
    scale_color_manual(values=cc) + 
    guides(colour = guide_legend(override.aes = list(size=2))) 
  
  return(scat)
}

S.abs <- function(df,xcol,ycol,x.tic,y.tic,cc,xlab,ylab){
  scat <- ggplot(df, aes(x=abs(df[[xcol]]),y=df[[ycol]])) + 
    geom_point(alpha=1, aes(color=df$Hypothesis), size=0.5) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=14, color="black"),
      text = element_text(size=14),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid = element_blank(), # get rid of major grid
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim=y.tic,xlim = x.tic) +
    scale_y_continuous(breaks = y.tic) + 
    scale_x_continuous(breaks = x.tic) +
    labs(x=xlab,y=ylab,color="Hypothesis") + 
    scale_color_manual(values=cc) + 
    guides(colour = guide_legend(override.aes = list(size=2))) 
  
  return(scat)
}


# col numbers for metadata 
hot <- 2 # alignment quality score 
seq <- 3 # total number of taxa
col <- 4 # sites in alignment
dist <- 5 # distinct patterns
pi <- 6 # PI sites
sing <- 7 # singleton sites
cons <- 8 # constant sites
chi2 <- 9 # number of sequences that failed composition chi2 test 
gaps <- 10 # number of sequences that contain 50% or more gaps/ambiguity
gl <- 11
hyp <- 12
bf <- 13 # max 180 for tox, all but 3 btw -50 and 20, 70 for tvs


yV <- gaps
yL <- "gaps"
#xV <- bf
#xL <- "ln(BF)"
#lines <- geom_vline(xintercept=c(-10,10),color=c("black"), linetype="dashed", size=0.5)
xV <- gl
xL <- "dGLS"
lines <- geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)


df <- AI.a
# Get max min for graph x = dgls y=bf
max(abs(df[[bf]]),na.rm = T)
max(abs(df[[gaps]]),na.rm = T)
x.lim <- 180 
x.t <- seq(-x.lim,x.lim,40)
y.lim <- 8
y.t <- seq(0,8,2)


#quartz()
# Set colors 
cc <- c('#1e7b59','#46d29f') # 3 up and 3 down from original
# 11 = gl, 13 = bf
  
AI.s <- S(AI.a,xV,yV,x.t,y.t,cc,xL,yL) +  lines

AI.s 

# Set colors 
cc <- c('#243a5b','#4975b6')

SA.s <- S(SA.a,xV,yV,x.t,y.t,cc,xL,yL) +  lines

SA.s

# Set colors 
cc <- c('#4d4d00','#cccc00')

SI.s <- S(SI.a,xV,yV,x.t,y.t,cc,xL,yL) +  lines


SI.s

cc = color_TP
TS.s <- S(TS.a,xV,yV,x.t,y.t,cc,xL,yL) +  lines

TS.s

SP <- ggarrange(TS.s, AI.s, SA.s, SI.s, ncol=1, nrow=4, align="h")
SP
#yL <- "gaps"
ggsave(paste(dataset,"_scatter_",xL,"_",yL,".pdf",sep=""), plot=SP,width = 6, height = 9, units = "in", device = 'pdf',bg = "transparent")




# Get max min for graph 

max(abs(min(mL$BF)),abs(max(mL$BF)),abs(min(mL$dGLS)),abs(max(mL$dGLS)))
limit <- 90
tic <- seq(-limit,limit,10)

# Names "TvS", "AIvSA", "AIvSI", "SAvAI", "SAvSI", "SIvAI", "SIvSA", "TvS_support"

color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "#C2DF23FF"


quartz()
# Set colors 
color_h0 <- color_AI
color_h0 <- color_TP

graph_general <- ggplot(mL, aes(x=BF,y=MEAN_COL_SCORE)) + 
  geom_point(alpha=0.5, color=color_h0, size=1) + theme_bw() + theme(panel.border = element_blank()) +
  theme_classic() + 
  theme(
    axis.text = element_text(size=16, color="black"),
    text = element_text(size=20),
    legend.position = "none", 
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(), # get rid of major grid
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim=tic,xlim = tic) +
  scale_y_continuous(breaks = tic) + 
  scale_x_continuous(breaks = tic) +
  labs(x='ln(BF)',y='dGLS')

graph_general

graph_custom <- graph_general + 
  geom_vline(xintercept=c(10,-10),color=c("gray"), linetype="dotted", size=0.5) +
  geom_hline(yintercept=c(0.5,-0.5),color=c("gray"), linetype="dotted", size=0.5) +
  ggtitle(paste(dataset,"\n",t,sep=""))
  #geom_abline(color=c("gray"), size=0.2, linetype="solid") + 
  

graph_custom

ggsave(paste(dataset,"_scatter_",hypoth,".pdf",sep=""), plot=graph_custom,width = 9.5, height = 6, units = "in", device = 'pdf',bg = "transparent")
