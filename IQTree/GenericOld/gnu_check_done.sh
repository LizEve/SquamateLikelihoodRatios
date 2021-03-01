#!/bin/bash

input='input_all.txt'
donelist='donelist_all.txt'

for f in *.nex
do 
run=`echo ${f%%.*}`
outfile1=`echo $run'/'$run'.treefile'`
outfile2=`echo $run'/'$run'.iqtree'`
# if output file doesnt exist, write to outfile 
if [[ -f $outfile1 && -f $outfile2 ]]
then
    echo $WDIR'/'$run'/' >> $donelist
else
    echo $WDIR'/'$run'/' >> $input
fi
done


'''
for f in *.nex
do 
run=`echo ${f%%.*}`
outfile=`echo $run'/'$run'.log'`
flen=`cat $outfile | wc -l`
echo $run
echo $flen
done
'''
