#rm(list=ls())
rm(z,zn,zp,mLBF,mLGL,gl,gl.s,G,g,c,bf.s,bf,b,B,anz,a,all.e,nLoci,dataset)

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stats)
library(tibble)
library(reshape)
# set dataset 

dataset <- "Reeder"

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

# Number of loci
nLoci <- length(mLGL$Locus)


# Transform datasets - stack two comparisons, and add column for bf and gl 
# Each dataset one hypothesis 

B <- mLBF %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType"))


G <- mLGL %>% 
  select(c("Locus","supportType","AIvSA","AIvSI","SAvSI","TvS")) %>% 
  melt(id=c("Locus","supportType")) 

# Smoosh together again
all.e <- left_join(G,B,by=c("Locus","variable"), suffix = c(".g", ".b")) 
a <- all.e %>% dplyr::rename(BF = value.b,GL = value.g) %>% 
  select(-c("supportType.b","supportType.g")) 

# Mark which are positive, so T = pos, F = neg, mark pos and neg significant, and zeros
#anz <- a %>% 
#  mutate(BF.p = BF > 0.00, GL.p = GL > 0.00,
#         BF.n = BF < 0.00, GL.n = GL < 0.00,
#         BF.sp = BF > 5, GL.sp = GL > 0.5,
#         BF.sn = BF < -5, GL.sn = GL < -0.5,
#         BF.z = BF == 0.00, GL.z = GL == 0.00) 

# Testing at 2 
anz <- a %>% 
  mutate(BF.p = BF > 0.00, GL.p = GL > 0.00,
         BF.n = BF < 0.00, GL.n = GL < 0.00,
         BF.sp = BF > 5, GL.sp = GL > 2,
         BF.sn = BF < -5, GL.sn = GL < -2,
         BF.z = BF == 0.00, GL.z = GL == 0.00) 

# Mean Mode Variance-----------------------------------------------------------------------------------------------------------------------------------------------------

zn <- anz %>% group_by(variable) %>%
  dplyr::summarise(BF.pos = sum(BF.p == TRUE),BF.neg = sum(BF.n == TRUE),BF.zero = sum(BF == 0.00),
                   GL.pos = sum(GL.p == TRUE),GL.neg = sum(GL.n == TRUE),GL.zero = sum(GL == 0.00)) 


# Add percents
zp <- zn %>% mutate_if(is.numeric, ~((./nLoci)*100)) %>% mutate_if(is.numeric,round,2)

z <- left_join(zn,zp,by=c("variable"),suffix = c(".n",".p"))

# Get stats for all positive and negative values
bf <- anz %>% 
  group_by(variable,BF.p,BF.n) %>% 
  dplyr::summarise(Total.Loci=n(),Mean=mean(BF),
             Median=median(BF), 
             Std=sd(BF),
             Max=max(BF),
             Min=min(BF)) %>% 
  mutate(metric = case_when(BF.p == TRUE ~ "BF, Pos",
                        BF.n == TRUE ~ "BF, Neg",
                        BF.n == FALSE & BF.p == FALSE ~ "BF, Zero"))

# Get stats for only significant numbers 
bf.s <- anz %>% 
  group_by(variable,BF.sp,BF.sn) %>% 
  dplyr::summarise(Total.Loci=n(),Mean=mean(BF),
            Median=median(BF), 
            Std=sd(BF),
            Max=max(BF),
            Min=min(BF)) %>% 
  mutate(metric = case_when(BF.sp == TRUE ~ "BF, Pos Sig",
                        BF.sn == TRUE ~ "BF, Neg Sig",
                        BF.sp==FALSE & BF.sn == FALSE ~ "BF, Not Sig"))
  

gl <- anz %>% 
  group_by(variable,GL.p,GL.n) %>% 
  dplyr::summarise(Total.Loci=n(),Mean=mean(GL),
            Median=median(GL), 
            Std=sd(GL),
            Max=max(GL),
            Min=min(GL))  %>% 
  mutate(metric = case_when(GL.p == TRUE ~ "GL, Pos",
                            GL.n == TRUE ~ "GL, Neg",
                            GL.n == FALSE & GL.p == FALSE ~ "GL, Zero"))

gl.s <- anz %>% 
  group_by(variable,GL.sp,GL.sn) %>% 
  dplyr::summarise(Total.Loci=n(),Mean=mean(GL),
            Median=median(GL), 
            Std=sd(GL),
            Max=max(GL),
            Min=min(GL)) %>% 
  mutate(metric = case_when(GL.sp == TRUE ~ "GL, Pos Sig",
                            GL.sn == TRUE ~ "GL, Neg Sig",
                            GL.sp == FALSE & GL.sn == FALSE ~ "GL, Not Sig"))
  
b <- bind_rows(bf, bf.s) %>% ungroup() %>%select(-c(BF.p, BF.n, BF.sp, BF.sn)) 
g <- bind_rows(gl, gl.s) %>% ungroup() %>%select(-c(GL.p, GL.n, GL.sp, GL.sn))

c <- bind_rows(b,g) %>% mutate(info=metric) %>% 
  tidyr::separate(metric,c("metric","value"),sep=", ") %>% 
  mutate(Percent.Loci =((Total.Loci/nLoci)*100)) %>% 
  mutate_if(is.numeric,round,2) %>%
  select(variable,metric,value,Percent.Loci,everything())


# reassign to dataset specific variable

assign(paste(dataset,"z",sep="."), z)
assign(paste(dataset,"c",sep="."), c)

#----Combine dataset-----------------------------------------------------------------------------------------------------------------------------------------------------


rm(z,zn,zp,mLBF,mLGL,gl,gl.s,G,g,c,bf.s,bf,b,B,anz,a,all.e,nLoci,dataset)



Burbrink.c <- Burbrink.c %>% add_column(Data.Set = "Burbrink", .before="variable")
Burbrink.z <- Burbrink.z %>% add_column(Data.Set = "Burbrink", .before="variable")

Singhal.c <- Singhal.c %>% add_column(Data.Set = "Singhal", .before="variable")
Singhal.z <- Singhal.z %>% add_column(Data.Set = "Singhal", .before="variable")

Reeder.c <- Reeder.c %>% add_column(Data.Set = "Reeder", .before="variable")
Reeder.z <- Reeder.z %>% add_column(Data.Set = "Reeder", .before="variable")

Streicher.c <- Streicher.c %>% add_column(Data.Set = "Streicher", .before="variable")
Streicher.z <- Streicher.z %>% add_column(Data.Set = "Burbrink", .before="variable")


c2 <- bind_rows(Burbrink.c,Streicher.c,Reeder.c,Singhal.c)
z2 <- bind_rows(Burbrink.z,Streicher.z,Reeder.z,Singhal.z)

#setwd("/Users/ChatNoir/Projects/Squam/scripts/Graphing/DataFiles/")
#save(c2,z2, file=paste("Calcs_countLociALL2.RData",sep=""))
#fname <- paste("/Users/ChatNoir/Projects/Squam/Graphs/ALL_hypothSupport2.xlsx",sep='')
#write.xlsx2(as.data.frame(c2), file=fname, sheetName="ALL2", row.names=FALSE)
#write.xlsx2(as.data.frame(z2), file=fname, sheetName="PositiveNegativeCounts2", append=TRUE, row.names=FALSE)


#c %>% group_by(variable) %>% dplyr::summarise(count=n())

#rm(z,zn,zp,mLBF,mLGL,gl,gl.s,G,g,c,bf.s,bf,b,B,anz,a,all.e)




##Depreciated -----------------------------------------------------------------------------------------------------------------------------------------------------
# 
# 
# #no zeros, plus metric column
# anzm <- a %>% filter(BF != 0, GL != 0) %>%
#   mutate(BF.direction = case_when(BF > 0 ~ "BF, Pos",
#                                   BF < 0 ~ "BF, Neg"),
#          BF.significance = case_when(BF >= 5 ~ "BF, Pos Sig",
#                                      BF <= -5 ~ "BF, Neg Sig")) %>%
#   mutate(GL.direction = case_when(GL > 0 ~ "GL, Pos",
#                                   GL < 0 ~ "GL, Neg"),
#          GL.significance = case_when(GL >= 0.5 ~ "GL, Pos Sig",
#                                      GL <= -0.5 ~ "GL, Neg Sig"))
# 
# direction <- anzm %>% group_by(variable,BF.direction,GL.direction) %>%
#   summarise(Total.Loci=n(),Mean.B=mean(BF), Median.B=median(BF), 
#             Std.B=sd(BF),Max.B=max(BF), Min.B=min(BF),Mean.G=mean(GL),
#             Median.G=median(GL), 
#             Std.G=sd(GL),
#             Max.G=max(GL),
#             Min.G=min(GL))
# 
# significant <- anzm %>% group_by(BF.significance,GL.significance) %>%
#   summarise(Total.Loci=n(),Mean.B=mean(BF), Median.B=median(BF), 
#             Std.B=sd(BF),Max.B=max(BF), Min.B=min(BF),Mean.G=mean(GL),
#             Median.G=median(GL), 
#             Std.G=sd(GL),
#             Max.G=max(GL),
#             Min.G=min(GL))