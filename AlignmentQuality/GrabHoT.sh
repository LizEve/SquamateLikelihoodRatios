#!/bin/bash
dataSet='Singhal'


working='/home/gmount/HoT/'$dataSet
storage='/media/ExtraDrive3/HoT/'$dataSet'_HoT'

# Run data gathering script 

# Output file 
outsummary=$working'/'$dataSet'_HoT.out'


# Gather info into data sheet 

echo -ne "Locus\t" > $outsummary
echo "MEAN_COL_SCORE" >> $outsummary

for f in $storage/*fasta
do
i=$(basename "$f" .fasta)
oldFile=$storage'/'$i'_cos_MFT/heads_HoT_msa.scr'
newFile=$working'/'$i'_heads_HoT_msa.scr'
cp $oldFile $newFile
echo -ne $i'\t' 		>> $outsummary
echo -ne `grep --text 'MEAN_COL_SCORE' $newFile | awk '{print $4}'`'\t'         >> $outsummary
echo "" >> $outsummary
done 