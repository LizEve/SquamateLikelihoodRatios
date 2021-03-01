#!/bin/bash
locusFasta=$1

working=`pwd`
storage='/media/ExtraDrive3/HoT/Singhal_HoT'

# get just the locus without extension
i="${locusFasta%.*}"

echo $i
pwd

# Set path variables - works on both mike and gorilla, same paths
COSPATH='/home/gmount/COS_LRM_v2.05/COS.pl'
MSASETPATH='/home/gmount/COS_LRM_v2.05/msa_set_score_v2.02/msa_set_score'

outfile=$storage"/"${i}_cos_MFT/heads_HoT_msa.scr
echo $outfile

if [ -f $outfile ] # if outfile exists
then
    echo $outfile" done"
else
    rsync -au $storage'/'$locusFasta $working
    wait
    # Run HoT scripts
    $COSPATH $i MFT nt $locusFasta `pwd`
    wait
    echo "step 1 done"
    $MSASETPATH $i"_cos_MFT/hot_H.fasta" $i"_cos_MFT/heads_HoT" -m $i"_cos_MFT/hot_T.fasta"
    wait
    echo "step 2 done"
    #Move back 
    rsync -au $i"_cos_MFT" $storage
    wait 
    rm -rf $i
    echo "moved and removed"
fi

