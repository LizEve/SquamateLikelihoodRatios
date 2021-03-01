#! /usr/bin/env python
import glob 
import csv 


# This script needs to be run in a folder with all *.out consel files, one for each locus. 
# Takes outfiles and concatenates them into one larger file, and adds a column for locus name. 
# Remove all other files that might have .out extentions, otherwise this will not work. 

# Name outfile and header
outFileName='Tox.csv' 
header=["LOCUS", "rank", "item", "obs", "au", "np", "bp", "pp", "kh", "sh", "wkh", "wsh"]
# Open outfile and write header 
newFile=open(outFileName, "w+")
wr = csv.writer(newFile)
wr.writerow(header)

# Suffix for outfiles from Consel 
suffix="tox3.out"

# Iterate through all outfiles in directory 
for f in glob.glob('./*%s' % suffix):
    # Get locus name 
    locus=f.strip(".").strip("/").split(".")[0].split("_")[0]
    # Open in file
    with open(f) as infile:
        # Read in all lines
        lines=infile.readlines()
        # Iterate through lines with data 
        for h in range(3,6):
            # Initiate a clean list with locus names
            cleanList=[locus]
            # Strip newline chars from line
            l=lines[h].strip()
            # Make line into list 
            dirtyList=list(filter(None,l.strip("#").split("  ")))
            # Clean extra spaces and characters from list
            cleanList.extend([y.strip(" ") for y in [x.strip("|") for x in dirtyList]])
            # Write list to master outfile 
            wr.writerow(cleanList)
    infile.close()