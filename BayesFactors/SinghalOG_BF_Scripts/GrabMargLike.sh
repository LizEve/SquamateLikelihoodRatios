#!/bin/bash
outsummary=$1
folderlist=$2
# Creates header line to specify marginal likelihoods collected
echo -ne "Locus\t" > $outsummary
echo -ne "Sclero\t" >> $outsummary
echo "ToxPoly\t" >> $outsummary



# While reading list of folders
while read f; do

# Grab locus name without the "/"
i=$(basename "$f" /)
# Add column with locus name, then marg LH values for each constraint
echo -ne $i'\t' 		>> $outsummary
echo -ne `grep --text Mean\: $i/$i".sclero.log" | awk '{print $2}'`'\t' 		>> $outsummary
echo -ne `grep --text Mean\: $i/$i".toxpoly.log" | awk '{print $2}'`'\t' 		>> $outsummary
echo "" >> ../$outsummary
done < $folderlist