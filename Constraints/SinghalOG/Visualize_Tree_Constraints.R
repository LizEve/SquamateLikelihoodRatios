library(ape)

setwd("/Users/ChatNoir/Projects/Squam/scripts/Constraints/SinghalOG")
t <- read.tree("scleroglossa.constraint")
y <- ladderize(t, right = TRUE)
pdf(file = "Sclero.pdf", width = 2, height = 10)
plot(y,cex=0.5,no.margin=TRUE)
dev.off()

t <- read.tree("toxicoferaP.constraint")
y <- ladderize(t, right = TRUE)
pdf(file = "ToxPoly.pdf", width = 2, height = 10)
plot(y,cex=0.5,no.margin=TRUE)
dev.off()
