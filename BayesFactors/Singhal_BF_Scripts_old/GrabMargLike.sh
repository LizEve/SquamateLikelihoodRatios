#!/bin/bash
outsummary=$1
folderlist=$2
# Creates header line to specify marginal likelihoods collected
echo -ne "Locus\t" > $outsummary
echo -ne "Sclero\t" >> $outsummary
echo -ne "ToxPoly\t" >> $outsummary
echo -ne "ToxAngIg\t" >> $outsummary
echo -ne "ToxSnAng\t" >> $outsummary
echo "ToxSnIg" >> $outsummary


# While reading list of folders
while read f; do
# Move into folder
cd $f
# Grab locus name without the "/"
i=$(basename "$f" /)
# Add column with locus name, then marg LH values for each constraint
echo -ne $i'\t' 		>> ../$outsummary
echo -ne `grep --text Mean\: $i".sclero.log" | awk '{print $2}'`'\t' 		>> ../$outsummary
echo -ne `grep --text Mean\: $i".toxpoly.log" | awk '{print $2}'`'\t' 		>> ../$outsummary
echo -ne `grep --text Mean\: $i".toxai.log" | awk '{print $2}'`'\t' 		>> ../$outsummary
echo -ne `grep --text Mean\: $i".toxsa.log" | awk '{print $2}'`'\t' 		>> ../$outsummary
echo -ne `grep --text Mean\: $i".toxsi.log" | awk '{print $2}'`'\t' 		>> ../$outsummary
echo "" >> ../$outsummary
cd ../
done < $folderlist