#!/bin/bash

runName="PRG"
input=input_${runName}.txt
donelist=donelist_${runName}.txt

JOBS_PER_NODE=8
TASK="$WDIR/task_iq_pbs"$runName".sh"
FILES="$WDIR/$input"

WDIR=`pwd`
cd $WDIR
rm $input
rm $donelist

for f in *.nex
do
r=$(basename "$f" .ntg.nex)
prefix=${r}.${constraint%%.*}
outfile1=$r/${prefix}.treefile
outfile2=$r/${prefix}.iqtree
#echo $outfile1
#echo $outfile2
# if output file doesnt exist, write to outfile
if [[ -f $outfile1 && -f $outfile2 ]]
then
    echo ${WDIR}/${r}/ >> $donelist
else
    echo ${WDIR}/${r}/ >> $input
fi
done

# Set up parallel call
# set log file, pass workers, whatever nodefile is, and working directory
PARALLEL="parallel --joblog gnu_$runName.log -j $WPN --wd $WDIR"

# Now call parallel and pass files and task. 
$PARALLEL -a $FILES sh $TASK {} 
