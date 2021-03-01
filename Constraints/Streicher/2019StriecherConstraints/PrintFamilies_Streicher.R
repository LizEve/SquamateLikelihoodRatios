library(readxl)
library(dplyr)
library(rJava)
library(xlsx)

inFile <- "/Users/ChatNoir/Google\ Drive/LSU/Chapter\ Squamate/Notes/Organizing.xlsx"
setwd("/Users/ChatNoir/Projects/Squam/scripts/Constraints")

# Read in xlsx file and name sheets. I got this from online somewhere. 
sheets <- readxl::excel_sheets(inFile)
x <- lapply(sheets, function(X) readxl::read_excel(inFile, sheet = X))
if(!FALSE) x <- lapply(x, as.data.frame)
names(x) <- sheets

# Put each sheet in it's own dataframe
burbrink <- as.data.frame(x[2])
streicher <- as.data.frame(x[3])
reederwiens <- as.data.frame(x[4])
singhal <- as.data.frame(x[5])

# Idk what this is 
#masterdf <- data.frame(ConstraintName=character(s),TaxaList=character(s),stringsAsFactors = FALSE)

### List families ##############################################
fams <- unique(streicher$StreicherTaxa.Family) 
# function to grab all taxa with family x
byfams <- function(x){
  streicher %>% filter(StreicherTaxa.Family==x) %>% 
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
  family <- unique(a$StreicherTaxa.Family)
  # get taxa 
  taxa <- a$StreicherTaxa.Original_names_of_sequences
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
fams <- unique(streicher$StreicherTaxa.Clade) 
# function to grab all taxa with family x
byfams <- function(x){
  streicher %>% filter(StreicherTaxa.Clade==x) %>% 
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
  family <- unique(a$StreicherTaxa.Clade)
  # get taxa 
  taxa <- a$StreicherTaxa.Original_names_of_sequences
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
fams <- unique(streicher$StreicherTaxa.SnakeClade) 
# function to grab all taxa with family x
byfams <- function(x){
  streicher %>% filter(StreicherTaxa.SnakeClade==x) %>% 
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
  family <- unique(a$StreicherTaxa.SnakeClade)
  # get taxa 
  taxa <- a$StreicherTaxa.Original_names_of_sequences
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
fams <- unique(streicher$StreicherTaxa.Scleroglossa) 
# function to grab all taxa with family x
byfams <- function(x){
  streicher %>% filter(StreicherTaxa.Scleroglossa==x) %>% 
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
  family <- unique(a$StreicherTaxa.Scleroglossa)
  # get taxa 
  taxa <- a$StreicherTaxa.Original_names_of_sequences
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
fams <- unique(streicher$StreicherTaxa.ToxPoly) 
# function to grab all taxa with family x
byfams <- function(x){
  streicher %>% filter(StreicherTaxa.ToxPoly==x) %>% 
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
  family <- unique(a$StreicherTaxa.ToxPoly)
  # get taxa 
  taxa <- a$StreicherTaxa.Original_names_of_sequences
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

### List ToxAI #############################################
fams <- 0
fams <- unique(streicher$StreicherTaxa.ToxAI) 
# function to grab all taxa with family x
byfams <- function(x){
  streicher %>% filter(StreicherTaxa.ToxAI==x) %>% 
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
  family <- unique(a$StreicherTaxa.ToxAI)
  # get taxa 
  taxa <- a$StreicherTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]

# Rename and clear 
dfToxAI <- df
df <- 0

### List ToxSA #############################################
fams <- 0
fams <- unique(streicher$StreicherTaxa.ToxSA) 
# function to grab all taxa with family x
byfams <- function(x){
  streicher %>% filter(StreicherTaxa.ToxSA==x) %>% 
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
  family <- unique(a$StreicherTaxa.ToxSA)
  # get taxa 
  taxa <- a$StreicherTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]

# Rename and clear 
dfToxSA <- df
df <- 0

### List ToxSI #############################################
fams <- 0
fams <- unique(streicher$StreicherTaxa.ToxSI) 
# function to grab all taxa with family x
byfams <- function(x){
  streicher %>% filter(StreicherTaxa.ToxSI==x) %>% 
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
  family <- unique(a$StreicherTaxa.ToxSI)
  # get taxa 
  taxa <- a$StreicherTaxa.Original_names_of_sequences
  if (length(taxa) > 1) {
    # Add family name to dataframe
    df$ConstraintName[n] <- as.character(family)
    # add taxa to data frame 
    df$TaxaList[n] <- list(taxa)
  }
}

df <- df[!apply(df == "", 1, all),]

# Rename and clear 
dfToxSI <- df
df <- 0

# Add them all together
new <- rbind(rbind(rbind(rbind(rbind(rbind(rbind(dfClade,dfFam),dfSnek),dfScler),dfToxPoly),dfToxAI),dfToxSA),dfToxSI)
# remove duplicate rows 
newer <- new[!duplicated(new),]
newer <- newer[- grep("Other", newer$ConstraintName),]

# print to file for IQ, only list, not name 
filename="Streicher_all_IQ.txt"
write(toString(newer[1,]), file=filename,append=FALSE)
for (i in 2:length(newer[,1])) {
  write(toString(newer[i,]), file=filename,append=TRUE)
}

# print to file for MB 
filename="Streicher_all.txt"
write(toString(newer[1,]), file=filename,append=FALSE)
for (i in 2:length(newer[,1])) {
  write(toString(newer[i,]), file=filename,append=TRUE)
}



# Get list of constraint names 
rownames(newer) <- c()
z <- 0
z <- as.character(newer$ConstraintName, rownames=FALSE)
x <- z[!duplicated(z)]
write(x,file="Streicher_ConstraintNames.txt",append=FALSE)

