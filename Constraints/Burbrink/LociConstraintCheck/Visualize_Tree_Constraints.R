library(ape)

setwd("/Users/ChatNoir/Projects/Squam/scripts/Constraints/Burbrink/LociConstraintCheck")
locus<-"T212_L14"
#locus<-'T212_L100'
tree<-'.treefile'
#tree<-''
t <- read.tree(paste(locus,".scleroglossa.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".sclero.pdf",sep = ''), width = 4, height = 17)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()

t <- read.tree(paste(locus,".toxicoferaP.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".toxPoly.pdf",sep = ''), width = 4, height = 17)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()

t <- read.tree(paste(locus,".toxicoferaAI.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".toxAI.pdf",sep = ''), width = 4, height = 17)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()

t <- read.tree(paste(locus,".toxicoferaSA.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".toxSA.pdf",sep = ''), width = 4, height = 17)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()

t <- read.tree(paste(locus,".toxicoferaSI.constraint",tree,sep = ''))
y <- ladderize(t, right = TRUE)
pdf(file = paste(locus,".toxSI.pdf",sep = ''), width = 4, height = 17)
plot(y,cex=0.5,no.margin=TRUE,use.edge.length = FALSE)
dev.off()