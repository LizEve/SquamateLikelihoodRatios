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

masterdf <- data.frame(ConstraintName=character(s),TaxaList=character(s),stringsAsFactors = FALSE)

### List families ##############################################
fams <- unique(burbrink$BurbrinkTaxa.Family) 
# function to grab all taxa with family x
byfams <- function(x){
  burbrink %>% filter(BurbrinkTaxa.Family==x) %>% 
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
  family <- unique(a$BurbrinkTaxa.Family)
  # get taxa 
  taxa <- a$BurbrinkTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

# Remove empty rows 
df <- df[!apply(df == "", 1, all),]

# Write names of each clade and taxa in clade
filename="Burbrink_Families.txt"
write(toString(df[1,]), file=filename,append=FALSE)
for (i in 2:length(df[,1])) {
  write(toString(df[i,]), file=filename,append=TRUE)
}

# List names of each constraint/clade
rownames(df) <- c()
z <- 0
z <- as.character(df$ConstraintName, rownames=FALSE)
write(z,file="Burbrink_ConstraintNames.txt",append=FALSE)

dfFam <- df
### List Clades #############################################
fams <- 0
fams <- unique(burbrink$BurbrinkTaxa.Clade) 
# function to grab all taxa with family x
byfams <- function(x){
  burbrink %>% filter(BurbrinkTaxa.Clade==x) %>% 
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
  family <- unique(a$BurbrinkTaxa.Clade)
  # get taxa 
  taxa <- a$BurbrinkTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]


filename="Burbrink_Clades.txt"
write(toString(df[1,]), file=filename,append=FALSE)
for (i in 2:length(df[,1])) {
  write(toString(df[i,]), file=filename,append=TRUE)
}

rownames(df) <- c()
z <- 0
z <- as.character(df$ConstraintName, rownames=FALSE)
write(z,file="Burbrink_ConstraintNames.txt",append=TRUE)

### List SnakeClade #############################################
fams <- 0
fams <- unique(burbrink$BurbrinkTaxa.SnakeClade) 
# function to grab all taxa with family x
byfams <- function(x){
  burbrink %>% filter(BurbrinkTaxa.SnakeClade==x) %>% 
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
  family <- unique(a$BurbrinkTaxa.SnakeClade)
  # get taxa 
  taxa <- a$BurbrinkTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]


filename="Burbrink_SnakeClade.txt"
write(toString(df[1,]), file=filename,append=FALSE)
for (i in 2:length(df[,1])) {
  write(toString(df[i,]), file=filename,append=TRUE)
}

rownames(df) <- c()
z <- 0
z <- as.character(df$ConstraintName, rownames=FALSE)
write(z,file="Burbrink_ConstraintNames.txt",append=TRUE)

### List Scleroglossa #############################################
fams <- 0
fams <- unique(burbrink$BurbrinkTaxa.Scleroglossa) 
# function to grab all taxa with family x
byfams <- function(x){
  burbrink %>% filter(BurbrinkTaxa.Scleroglossa==x) %>% 
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
  family <- unique(a$BurbrinkTaxa.Scleroglossa)
  # get taxa 
  taxa <- a$BurbrinkTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]


filename="Burbrink_Scleroglossa.txt"
write(toString(df[1,]), file=filename,append=FALSE)
for (i in 2:length(df[,1])) {
  write(toString(df[i,]), file=filename,append=TRUE)
}

rownames(df) <- c()
z <- 0
z <- as.character(df$ConstraintName, rownames=FALSE)
write(z,file="Burbrink_ConstraintNames.txt",append=TRUE)

### List ToxPoly #############################################
fams <- 0
fams <- unique(burbrink$BurbrinkTaxa.ToxPoly) 
# function to grab all taxa with family x
byfams <- function(x){
  burbrink %>% filter(BurbrinkTaxa.ToxPoly==x) %>% 
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
  family <- unique(a$BurbrinkTaxa.ToxPoly)
  # get taxa 
  taxa <- a$BurbrinkTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]


filename="Burbrink_ToxPoly.txt"
write(toString(df[1,]), file=filename,append=FALSE)
for (i in 2:length(df[,1])) {
  write(toString(df[i,]), file=filename,append=TRUE)
}

rownames(df) <- c()
z <- 0
z <- as.character(df$ConstraintName, rownames=FALSE)
write(z,file="Burbrink_ConstraintNames.txt",append=TRUE)

### List ToxAI #############################################
fams <- 0
fams <- unique(burbrink$BurbrinkTaxa.ToxAI) 
# function to grab all taxa with family x
byfams <- function(x){
  burbrink %>% filter(BurbrinkTaxa.ToxAI==x) %>% 
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
  family <- unique(a$BurbrinkTaxa.ToxAI)
  # get taxa 
  taxa <- a$BurbrinkTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]


filename="Burbrink_ToxAI.txt"
write(toString(df[1,]), file=filename,append=FALSE)
for (i in 2:length(df[,1])) {
  write(toString(df[i,]), file=filename,append=TRUE)
}

rownames(df) <- c()
z <- 0
z <- as.character(df$ConstraintName, rownames=FALSE)
write(z,file="Burbrink_ConstraintNames.txt",append=TRUE)

### List ToxSA #############################################
fams <- 0
fams <- unique(burbrink$BurbrinkTaxa.ToxSA) 
# function to grab all taxa with family x
byfams <- function(x){
  burbrink %>% filter(BurbrinkTaxa.ToxSA==x) %>% 
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
  family <- unique(a$BurbrinkTaxa.ToxSA)
  # get taxa 
  taxa <- a$BurbrinkTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]


filename="Burbrink_ToxSA.txt"
write(toString(df[1,]), file=filename,append=FALSE)
for (i in 2:length(df[,1])) {
  write(toString(df[i,]), file=filename,append=TRUE)
}


rownames(df) <- c()
z <- 0
z <- as.character(df$ConstraintName, rownames=FALSE)
write(z,file="Burbrink_ConstraintNames.txt",append=TRUE)


### List ToxSI #############################################
fams <- 0
fams <- unique(burbrink$BurbrinkTaxa.ToxSI) 
# function to grab all taxa with family x
byfams <- function(x){
  burbrink %>% filter(BurbrinkTaxa.ToxSI==x) %>% 
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
  family <- unique(a$BurbrinkTaxa.ToxSI)
  # get taxa 
  taxa <- a$BurbrinkTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]


filename="Burbrink_ToxSI.txt"
write(toString(df[1,]), file=filename,append=FALSE)
for (i in 2:length(df[,1])) {
  write(toString(df[i,]), file=filename,append=TRUE)
}

rownames(df) <- c()
z <- 0
z <- as.character(df$ConstraintName, rownames=FALSE)
write(z,file="Burbrink_ConstraintNames.txt",append=TRUE)
