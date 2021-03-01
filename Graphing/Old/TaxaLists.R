rm(list=ls())
library(dplyr)
library(tidyr)
library(ape)

# set dataset 
# read.nexus.data(paste("/Users/ChatNoir/Projects/Squam/Singhal/SinghalOG_Data/loci_files/",l,".nex",sep="")) 

dataset <- "Singhal"

setwd(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset, sep=''))

# Read in data
load(paste("/Users/ChatNoir/Projects/Squam/Graphs/",dataset,"/Calcs_",dataset,".RData", sep=''))

# Change from 2ln(BF) to ln(BF) by dividing all comparisons by 2 7-10 - not including average calcs
names(mLBF[7:10])
mLBF[7:10] <- mLBF[7:10]/2
# for singhal OG 
#names(mLBF[4])
#mLBF[4] <- mLBF[4]/2
metaData <- c("MEAN_COL_SCORE","Sequences","Columns","Dist_Pat","Pars_Info","Sing_Sites","Cons_Sites" ,"Chi2_Fail","Gaps_Ambig")
G <- mLGL %>% select(c(Locus,TvS))
B <- mLBF %>% select(c(Locus,TvS))
all <- left_join(B, G, by=c("Locus"), suffix = c(".g", ".b"))
rownames(all) <- all[,1]
allbkp <- all
# Add genome taxa columns 
genomes <- c("Alligator_mississippiensis", "Chrysemys_picta", "Eublepharis_macularius", "Gekko_japonicus", "Gopherus_evgoodei", "Lacerta_bilineata", "Lacerta_viridis", "Paroedura_picta", "Podarcis_muralis", "Salvator_meriana", "Sphenodon_punctatus1", "Sphenodon_punctatus2", "Zootoca_vivipara", "Anolis_carolinensis", "Crotalus_horridus", "Deinagkistrodon_actutus", "Dopasia_gracilis", "Hydrophis_melanocephalus", "Laticauda_colubrina", "Ophiophagus_hannah", "Pogona_vitticeps", "Protobothrops_mucrosquamatus", "Python_bivittatus", "Shinisaurus_crocodilurus", "Varanus_komodoensis")
for (i in genomes)
  all[,i] <- NA

# Add taxa from genome 

for (l in all$Locus){
  print(l)
  # Read in nexus file
  locus <- read.nexus.data(paste("/Users/ChatNoir/Projects/Squam/Singhal/loci_files/",l,".nex",sep="")) 
  # Get taxa names
  taxa <- names(locus)
  # Add a 1 for each taxa for this locus 
  for (c in taxa){
    print(c)
    all[l,c] <- 1
    print(all[l,c])
  }
}

# save cuz that took a few minutes
save(all,file=paste("GenomeTaxaList_",dataset,".RData",sep=""))

# Replace Na with 0
all[is.na(all)] = 0

# Remove columns accidentally added 
xlall.bkp <- all
all <- all.bkp

# How many genome taxa for each locus 
all <- all %>% mutate(num.gen = rowSums(.[4:28]))

# graph scatter 

## Scatter graph, colored by number of genome taxa
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
S <- function(df,xcol,ycol,x.tic,y.tic,cc,clab,xlab,ylab){
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]])) + 
    geom_point(alpha=1, aes(color=cc), size=0.95) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=12, color="black"),
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
    labs(x=xlab,y=ylab,color=clab)
    #scale_colour_gradient2()
    #scale_color_manual(values=cc) + 
    #guides(colour = guide_legend(override.aes = list(size=2))) 
  
  return(scat)
}


# Compare bf and gl vs genome taxa
cc <- all$num.gen
y.t <- seq(0,25,5)
x.t <- seq(-75,75,20)
a <- S(all,29,2,y.t,x.t,cc,'genome taxa','genome taxa','dGLS') + 
  geom_hline(yintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_viridis_c(option="plasma")
a
x.t <- seq(-85,85,20)
b <- S(all,29,3,y.t,x.t,cc,'genome taxa','genome taxa', 'ln(BF)') + 
  geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) + 
  scale_color_viridis_c(option="plasma")
b
ab <- ggarrange(a,b, ncol=2, nrow=1, align="v")
#ggsave("Singhal_dglandbfvgenometaxa_scatter.pdf", plot=ab,width = 9, height = 6, units = "in", device = 'pdf',bg = "transparent")

# bf v dgls colored by number of taxa

df <- all
# Get max min for graph x = dgls y=bf
max(abs(df[[2]]))
x.lim <- 75 # 50, 45
x.t <- seq(-x.lim,x.lim,20)
max(abs(df[[3]]))
y.lim <- 85 #1000, 120
y.t <- seq(-y.lim,y.lim,20)

# Set colors 
cc <- c('#1e7b59','#46d29f') # 3 up and 3 down from original

a <- S(all,2,3,x.t,y.t,cc,'dGLS','ln(BF)') + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)

#quartz()
a <- a + scale_color_viridis_c(option="plasma")
a
#ggsave("Singhal_dglsvbf_genometaxa_scatter.pdf", plot=a,width = 7, height = 6, units = "in", device = 'pdf',bg = "transparent")

# bf v dgls colored by presence/absense of taxa 

#a1 <- S(all,2,3,x.t,y.t,all$Sphenodon_punctatus1,"spheno1",'dGLS','ln(BF)') + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
#  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)
#a1

# Change binary present/absent to words and factor. in new data frame
allf <- all
for (x in seq(4,28,by=1)){
  allf[,x] <- gsub('1','present',allf[,x])
  allf[,x] <- gsub('0','absent',allf[,x])
  allf[,x] <- factor(allf[,x])
}


plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=df[,label_column])) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black",bins=100) +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=names(df)[label_column]))
}

#quartz()
out.cols <- c(4,5,8,14,15)
ig.cols <- c(17,24)
ang.cols <- c(20,27,28)
sna.cols <- c(25,26,18,19,21,22,23)
other.cols <- c(6,7,9,10,11,12,13,16)

make_hist = function (x){
  plot_multi_histogram(allf, 'TvS.b',x) 
}

myplots1 <- lapply(ig.cols,make_hist)
myplots2 <- lapply(ang.cols,make_hist)
myplots <- c(myplots1,myplots2)
myplots <- lapply(other.cols,make_hist)
ab <- ggarrange(plotlist=myplots, ncol=1, nrow=length(myplots), align="v")
ab
ggsave("Singhal_bf_other_hist.pdf", plot=ab,width = 6, height = 9, units = "in", device = 'pdf',bg = "transparent")

# Looks like sphenodon2 and varanus are crappo new scatter 

S <- function(df,xcol,ycol,x.tic,y.tic,cc,clab,xlab,ylab){
  scat <- ggplot(df, aes(x=df[[xcol]],y=df[[ycol]])) + 
    geom_point(alpha=0.5, aes(color=cc), size=1.5) + theme_bw() + theme(panel.border = element_blank()) +
    theme_classic() + 
    theme(
      axis.text = element_text(size=12, color="black"),
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
    labs(x=xlab,y=ylab,color=clab)
  return(scat)
}


# Make column to color loci with v,s,vs,n
temp <- allf
bt <- allf %>% mutate(bad.tax = case_when(Sphenodon_punctatus2 == 'present' & Varanus_komodoensis == 'present' ~ '1. both',
                                          Sphenodon_punctatus2 == 'present' & Varanus_komodoensis == 'absent' ~ '2. spheno2',
                                          Varanus_komodoensis == 'present' & Sphenodon_punctatus2 == 'absent' ~ '3. varanus')) 

bt$bad.tax[is.na(bt$bad.tax)] <- "4. neither"
h1 <- plot_multi_histogram(bt, 'TvS.g',31) 
h1
# Scatter plot colored by taxa presence
# Get max min for graph x = dgls y=bf
x.t <- seq(-75,75,20)
y.t <- seq(-85,85,20)

# Set colors 
cc <- bt$bad.tax
a1 <- S(bt,2,3,x.t,y.t,cc,'bad.tax','dGLS','ln(BF)') + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)
a1 <- a1 + scale_color_manual(values=c("green4","blue4","purple4",'gray81'))
a1 <- a1 + scale_color_manual(values=c("green4","gray81","gray81",'gray81'))
a1
# try with random taxa 

rt <- allf %>% mutate(bad.tax = case_when(Chrysemys_picta == 'present' & Dopasia_gracilis == 'present' ~ '1. both',
                                          Chrysemys_picta == 'present' ~ '2. chysemys',
                                          Dopasia_gracilis == 'present' ~ '3. dopasia')) 

rt$bad.tax[is.na(rt$bad.tax)] <- "4. neither"

# Set colors 
cc <- rt$bad.tax
a2 <- S(rt,2,3,x.t,y.t,cc,'bad.tax','dGLS','ln(BF)') + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)
a2 <- a2 + scale_color_manual(values=c("blue4",'gray81',"purple4","green4"))
a2 <- a2 + scale_color_manual(values=c("gray81","gray81",'gray81',"green4"))
h2 <- plot_multi_histogram(rt, 'TvS.g',31) 
h2
# Another time with random taxa
rt <- allf %>% mutate(bad.tax = case_when(Sphenodon_punctatus1 == 'present' & Shinisaurus_crocodilurus == 'present' ~ '1. both',
                                          Sphenodon_punctatus1 == 'present' ~ '2. spheno1',
                                          Shinisaurus_crocodilurus == 'present' ~ '3. shinis')) 

rt$bad.tax[is.na(rt$bad.tax)] <- "4. neither"

# Set colors 
cc <- rt$bad.tax
a3 <- S(rt,2,3,x.t,y.t,cc,'bad.tax','dGLS','ln(BF)') + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)
a3 <- a3 + scale_color_manual(values=c("green4","blue4",'gray81',"purple4"))
a3 <- a3 + scale_color_manual(values=c("green4","gray81",'gray81',"gray81"))
a3
a2

h3 <- plot_multi_histogram(rt, 'TvS.g',31) 
h3
ab <- ggarrange(a1,a2,a3, ncol=1, nrow=3, align="v")
ab

h <- ggarrange(h1,h2,h3, ncol=1, nrow=3, align="v")
h
#ggsave("Singhal_badtaxa_DGLS_scatter.pdf", plot=h,width = 6, height = 5, units = "in", device = 'pdf',bg = "transparent")


# only turtle as outgroup give different answers?

ot.all <- allf %>% 
  #filter(Varanus_komodoensis == 'absent') %>%
  #filter(Sphenodon_punctatus2 == 'absent') %>%
  filter(Sphenodon_punctatus2 == 'absent') %>%
  filter(Sphenodon_punctatus2 == 'absent') %>%
  mutate(ELT = case_when(Chrysemys_picta == 'present' & Gopherus_evgoodei == 'present' ~ 'both',
                                    Chrysemys_picta == 'present' ~ 'chrysemys',
                                    Gopherus_evgoodei == 'present' ~ 'gopherus')) 
ot$ELT[is.na(ot$ELT)] <- "neither"
ot.all$ELT[is.na(ot.all$ELT)] <- "neither"

cc <- ot$ELT
t3 <- S(ot,2,3,x.t,y.t,cc,'ELT','dGLS','ln(BF)') + geom_hline(yintercept=c(-5,5),color=c("black"), linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(-0.5,0.5),color=c("black"), linetype="dashed", size=0.5)
t3 <- t3 + scale_color_manual(values=c("green4","blue4","purple4",'yellow4'))

quartz()
plot_multi_histogram(ot, 'TvS.b',31) 
plot_multi_histogram(ot.all, 'TvS.b',31) 
