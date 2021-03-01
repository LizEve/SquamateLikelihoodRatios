library(readxl)
library(dplyr)
library(rJava)
library(xlsx)

inFile <- "/Users/ChatNoir/Google\ Drive/LSU/Chapter\ Squamate/Notes/Organizing.xlsx"
setwd("/Users/ChatNoir/Projects/Squam/scripts/Constraints/SinghalOG")

# Read in xlsx file and name sheets. I got this from online somewhere. 
sheets <- readxl::excel_sheets(inFile)
x <- lapply(sheets, function(X) readxl::read_excel(inFile, sheet = X))
if(!FALSE) x <- lapply(x, as.data.frame)
names(x) <- sheets

# Put each sheet in it's own dataframe
singhalog <- as.data.frame(x[7])

# Idk what this is 
#masterdf <- data.frame(ConstraintName=character(s),TaxaList=character(s),stringsAsFactors = FALSE)

### List families ##############################################
fams <- unique(singhalog$SinghalOGTaxa.Family) 
# function to grab all taxa with family x
byfams <- function(x){
  singhalog %>% filter(SinghalOGTaxa.Family==x) %>% 
    ungroup %>%
    as.data.frame()
}
# List of dataframes, parsed by family 
l <- lapply(fams,byfams)

# get number of families and create new dataframe
s <- length(fams) 
f <- seq.int(1,s,by=1) 
df <- data.frame(ConstraintName=character(s),TaxaList=character(s),stringsAsFactors = FALSE)

# iterate through all numbers in total number of families
for (n in f){
  # get dataframe of single family and taxa
  a <- l[[n]]
  # extract the family name
  family <- unique(a$SinghalOGTaxa.Family)
  # get taxa 
  taxa <- a$SinghalOGTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

# Remove empty rows 
df <- df[!apply(df == "", 1, all),]

# Rename and clear 
dfFam <- df
df <- 0

### List Clades #############################################
fams <- 0
fams <- unique(singhalog$SinghalOGTaxa.Clade) 
# function to grab all taxa with family x
byfams <- function(x){
  singhalog %>% filter(SinghalOGTaxa.Clade==x) %>% 
    ungroup %>%
    as.data.frame()
}
# List of dataframes, parsed by family 
l <- lapply(fams,byfams)

# get number of families and create new dataframe
s <- length(fams)
f <- seq(1,s,by=1)
df <- data.frame(ConstraintName=character(s),TaxaList=character(s),stringsAsFactors = FALSE)

# iterate through all numbers in total number of families
for (n in f){
  # get dataframe of single family and taxa
  a <- l[[n]]
  # extract the family name
  family <- unique(a$SinghalOGTaxa.Clade)
  # get taxa 
  taxa <- a$SinghalOGTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]

# Rename and clear 
dfClade <- df
df <- 0


### List SnakeClade #############################################
fams <- 0
fams <- unique(singhalog$SinghalOGTaxa.SnakeClade) 
# function to grab all taxa with family x
byfams <- function(x){
  singhalog %>% filter(SinghalOGTaxa.SnakeClade==x) %>% 
    ungroup %>%
    as.data.frame()
}
# List of dataframes, parsed by family 
l <- lapply(fams,byfams)

# get number of families and create new dataframe
s <- length(fams)
f <- seq(1,s,by=1)
df <- data.frame(ConstraintName=character(s),TaxaList=character(s),stringsAsFactors = FALSE)

# iterate through all numbers in total number of families
for (n in f){
  # get dataframe of single family and taxa
  a <- l[[n]]
  # extract the family name
  family <- unique(a$SinghalOGTaxa.SnakeClade)
  # get taxa 
  taxa <- a$SinghalOGTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]

# Rename and clear 
dfSnek <- df
df <- 0

### List Scleroglossa #############################################
fams <- 0
fams <- unique(singhalog$SinghalOGTaxa.Scleroglossa) 
# function to grab all taxa with family x
byfams <- function(x){
  singhalog %>% filter(SinghalOGTaxa.Scleroglossa==x) %>% 
    ungroup %>%
    as.data.frame()
}
# List of dataframes, parsed by family 
l <- lapply(fams,byfams)

# get number of families and create new dataframe
s <- length(fams)
f <- seq(1,s,by=1)
df <- data.frame(ConstraintName=character(s),TaxaList=character(s),stringsAsFactors = FALSE)

# iterate through all numbers in total number of families
for (n in f){
  # get dataframe of single family and taxa
  a <- l[[n]]
  # extract the family name
  family <- unique(a$SinghalOGTaxa.Scleroglossa)
  # get taxa 
  taxa <- a$SinghalOGTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]

# Rename and clear 
dfScler <- df
df <- 0


### List ToxPoly #############################################
fams <- 0
fams <- unique(singhalog$SinghalOGTaxa.ToxPoly) 
# function to grab all taxa with family x
byfams <- function(x){
  singhalog %>% filter(SinghalOGTaxa.ToxPoly==x) %>% 
    ungroup %>%
    as.data.frame()
}
# List of dataframes, parsed by family 
l <- lapply(fams,byfams)

# get number of families and create new dataframe
s <- length(fams)
f <- seq(1,s,by=1)
df <- data.frame(ConstraintName=character(s),TaxaList=character(s),stringsAsFactors = FALSE)

# iterate through all numbers in total number of families
for (n in f){
  # get dataframe of single family and taxa
  a <- l[[n]]
  # extract the family name
  family <- unique(a$SinghalOGTaxa.ToxPoly)
  # get taxa 
  taxa <- a$SinghalOGTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]
# Rename and clear 
dfToxPoly <- df
df <- 0

# Add them all together
new <- rbind(rbind(rbind(rbind(dfClade,dfFam),dfSnek),dfScler),dfToxPoly)
# remove duplicate rows 
newer <- new[!duplicated(new),]
newer <- newer[- grep("Other", newer$ConstraintName),]

# added for IQ constraints 
newer$TaxaListClean <- gsub(pattern = '"',replacement ='',newer$TaxaList,perl=TRUE)



######### Print constraint files 



# Make dataset for each constraint. List clades to be removed. 
Sclero <- newer %>% filter(!ConstraintName %in% c('ToxPoly'))
ToxPoly <- newer %>% filter(!ConstraintName %in% c('Sclero'))


# print to file for IQ, only list, not name 
makeFile = function(df){
  fname = paste('SinghalOG_',deparse(substitute(df)),'_iq.txt',sep='')
  write.table(t(df$TaxaListClean), col.names = FALSE, row.names = FALSE,file = fname)
  print(fname)
}

makeFile(Sclero)
makeFile(ToxPoly)



# print to file for MB 
filename="SinghalOG_all_bb.txt"
newbb <- new[!duplicated(new),]
newbb <- newbb[- grep("Other", newbb$ConstraintName),]

write(toString(newbb[1,]), file=filename,append=FALSE)
for (i in 2:length(newbb[,1])) {
  write(toString(newbb[i,]), file=filename,append=TRUE)
}



# Get list of constraint names 
rownames(newbb) <- c()
z <- 0
z <- as.character(newbb$ConstraintName, rownames=FALSE)
x <- z[!duplicated(z)]
write(x,file="SinghalOG_ConstraintNames.txt",append=FALSE)

