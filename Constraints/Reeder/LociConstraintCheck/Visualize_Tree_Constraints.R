library(ape)

setwd("/Users/ChatNoir/Projects/Squam/scripts/Constraints/Reeder/LociConstraintCheck")
locus<-"Reeder_DNA_BACH1"
#locus<-'Reeder_DNA_BHLHB2'
tree<-'.treefile'
t <- read.tree(paste(locus,".scleroglossa.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".sclero.pdf",sep = ''), width = 2, height = 10)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()

t <- read.tree(paste(locus,".toxicoferaP.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".toxPoly.pdf",sep = ''), width = 2, height = 10)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()

t <- read.tree(paste(locus,".toxicoferaAI.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".toxAI.pdf",sep = ''), width = 2, height = 10)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()

t <- read.tree(paste(locus,".toxicoferaSA.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".toxSA.pdf",sep = ''), width = 2, height = 10)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()

t <- read.tree(paste(locus,".toxicoferaSI.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".toxSI.pdf",sep = ''), width = 2, height = 10)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()