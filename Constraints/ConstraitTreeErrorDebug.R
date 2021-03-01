library(ape)

setwd("/Users/ChatNoir/Projects/Squam/Streicher/StreicherData/uce-497")
t <- read.tree("uce-497.scleroglossa.constraint")
y <- ladderize(t, right = TRUE)
#pdf(file = "uce-497.scleroglossa.constraint.pdf", width = 2, height = 10)
quartz()
plot(y,cex=0.5,no.margin=TRUE)
#dev.off()


t <- read.tree(text='(sphenodon_sp,(strophurus_sp,saltorius_sp,gekko_sp,coleonyx_variegatus,lialis_sp,(cordylosaurus_sp,cordylus_sp,cricosaura_sp,plestiodon_fasciatus,tiliqua_sp),(aspidoscelis_tigris,tupinambus_sp),(bipes_sp,rhineura_sp),(micrurus_fulvius,python_molurus_2),varanus_exanthematicus));')
y <- ladderize(t, right = TRUE)
quartz()
plot(y,cex=0.5,no.margin=TRUE)