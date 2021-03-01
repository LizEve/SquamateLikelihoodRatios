library(readxl)
library(dplyr)
library(rJava)
library(xlsx)

inFile <- "/Users/ChatNoir/Google\ Drive/LSU/Chapter\ Squamate/Notes/Organizing.xlsx"
setwd("/Users/ChatNoir/Projects/Squam/Constraints")

# Read in xlsx file and name sheets. I got this from online somewhere. 
sheets <- readxl::excel_sheets(inFile)
x <- lapply(sheets, function(X) readxl::read_excel(inFile, sheet = X))
if(!FALSE) x <- lapply(x, as.data.frame)
names(x) <- sheets

# Put each sheet in it's own dataframe
burbrink <- as.data.frame(x[2])
streicher <- as.data.frame(x[4])
reederwiens <- as.data.frame(x[5])
singhal <- as.data.frame(x[3])

# Count taxa per clade 

burbrink.clade.count <- burbrink %>% group_by(BurbrinkTaxa.Bigger.Clade,BurbrinkTaxa.Clade) %>% summarize(CladeCount.Burbrink=n())
streicher.clade.count <- streicher %>% group_by(StreicherTaxa.Bigger.Clade, StreicherTaxa.Clade) %>% summarize(CladeCount.Streicher=n())
reederwiens.clade.count <- reederwiens %>% group_by(ReederTaxa.Bigger.Clade,ReederTaxa.Clade) %>% summarize(CladeCount.Reeder=n())
singhal.clade.count <- singhal %>% group_by(SinghalTaxa.Bigger.Clade,SinghalTaxa.Clade) %>% summarize(CladeCount.Singhal=n())


merge(burbrink.clade.count,streicher.clade.count, 
      by.x=c("BurbrinkTaxa.Bigger.Clade","BurbrinkTaxa.Clade"), by.y=c("StreicherTaxa.Bigger.Clade","StreicherTaxa.Clade"), all = TRUE
)

# Join all together to get matrix of taxa per clade 
# burbrink has most clades, so use it as base 
clade.count <- merge(
  merge(
    merge(burbrink.clade.count,streicher.clade.count, 
          by.x=c("BurbrinkTaxa.Bigger.Clade","BurbrinkTaxa.Clade"), by.y=c("StreicherTaxa.Bigger.Clade","StreicherTaxa.Clade"), all = TRUE
          ), 
    reederwiens.clade.count, by.x=c("BurbrinkTaxa.Bigger.Clade","BurbrinkTaxa.Clade"), by.y=c("ReederTaxa.Bigger.Clade","ReederTaxa.Clade"), all = TRUE
    ), 
  singhal.clade.count, by.x=c("BurbrinkTaxa.Bigger.Clade","BurbrinkTaxa.Clade"), by.y=c("SinghalTaxa.Bigger.Clade","SinghalTaxa.Clade"), all = TRUE
)


# Count taxa per Family

burbrink.family.count <- burbrink %>% group_by(BurbrinkTaxa.Clade,BurbrinkTaxa.Family) %>% summarize(FamilyCount.Burbrink=n())
streicher.family.count <- streicher %>% group_by(StreicherTaxa.Clade,StreicherTaxa.Family) %>% summarize(FamilyCount.Streicher=n())
reederwiens.family.count <- reederwiens %>% group_by(ReederTaxa.Clade,ReederTaxa.Family) %>% summarize(FamilyCount.Reeder=n())
singhal.family.count <- singhal %>% group_by(SinghalTaxa.Clade,SinghalTaxa.Family) %>% summarize(FamilyCount.Singhal=n())

# Join all together to get matrix of taxa per family 
# burbrink has most familys, so use it as base 
family.count <- merge(
  merge(
    merge(burbrink.family.count,streicher.family.count, 
          by.x=c("BurbrinkTaxa.Clade","BurbrinkTaxa.Family"), by.y=c("StreicherTaxa.Clade","StreicherTaxa.Family"), all = TRUE
    ), 
    reederwiens.family.count, by.x=c("BurbrinkTaxa.Clade","BurbrinkTaxa.Family"), by.y=c("ReederTaxa.Clade","ReederTaxa.Family"), all = TRUE
  ), 
  singhal.family.count, by.x=c("BurbrinkTaxa.Clade","BurbrinkTaxa.Family"), by.y=c("SinghalTaxa.Clade","SinghalTaxa.Family"), all = TRUE
)


# Save as xlsx file 

write.xlsx(clade.count, file="CladeCounts.xlsx", sheetName="clade.count", row.names=FALSE)
write.xlsx(family.count, file="CladeCounts.xlsx", sheetName="family.count", append=TRUE, row.names=FALSE)


