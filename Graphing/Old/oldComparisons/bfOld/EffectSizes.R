library(ggplot2)
library(tidyverse)

setwd(getwd())
getwd()
uce314 <- "uce-314/results_GTRG/inference_pvalues_effectsizes_uce-314.csv"
uce12 <- "uce-12/results_GTRG/inference_pvalues_effectsizes_uce-12.csv"
ahe13 <- "AHE-L13/results_GTRG/inference_pvalues_effectsizes_AHE-L13.csv"

z <- read.csv(uce314, header=TRUE)
z$Locus <- rep("uce314", length(z$Statistic))
y <- read.csv(uce12, header=TRUE)
y$Locus <- rep("uce12", length(y$Statistic))
x <- read.csv(ahe13, header=TRUE)
x$Locus <- rep("ahe13", length(x$Statistic))

p <- subset(z, select = c("Statistic","Locus","Effect.Size"))
pp <- bind_rows(z,y)
ppes <- bind_rows(pp,x) 


a <- ggplot(ppes, aes(fill=Statistic, y=Effect.Size, x=Locus)) + 
  geom_bar(position="dodge", stat="identity")
ggsave("EffectSize_Example.pdf", plot=a, width = 15, height = 10, units = "cm")
