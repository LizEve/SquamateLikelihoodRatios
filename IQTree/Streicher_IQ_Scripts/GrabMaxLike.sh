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
echo -ne `grep --text "Log-likelihood of the tree\:" *scleroglossa.*iqtree | awk '{print $5}'`'\t'         >> ../$outsummary
echo -ne `grep --text "Log-likelihood of the tree\:" *toxicoferaP.*iqtree | awk '{print $5}'`'\t'         >> ../$outsummary
echo -ne `grep --text "Log-likelihood of the tree\:" *toxicoferaAI.*iqtree | awk '{print $5}'`'\t'         >> ../$outsummary
echo -ne `grep --text "Log-likelihood of the tree\:" *toxicoferaSA.*iqtree | awk '{print $5}'`'\t'         >> ../$outsummary
echo -ne `grep --text "Log-likelihood of the tree\:" *toxicoferaSI.*iqtree | awk '{print $5}'`'\t'         >> ../$outsummary
echo "" >> ../$outsummary
cd ../
done < $folderlist
