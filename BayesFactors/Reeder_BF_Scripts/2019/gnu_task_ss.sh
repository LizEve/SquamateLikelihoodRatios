#!/bin/bash
# Move into gene folder
cd $1
WDIR=`pwd`
run=`echo $(basename $WDIR )`

function runMB {
    run=$1
    outfile=$2
    bb=$3
    if [ -f $outfile ] # if outfile exists
    then
        if [[ $(wc -l < "$outfile") -eq 57 ]]
        then
            echo $outfile" done"
        else
            /home/gmount/MrBayes-3.2.7a/src/mb $bb # run mrbayes if output isnt long enough
        fi
    else
    /home/gmount/MrBayes-3.2.7a/src/mb $bb # run mrbayes if output doesnt exists. 
    fi
}


bb=`echo $run'_sclero.bb'`
outfile=${run}.sclero.ss
runMB $run $outfile $bb

bb=`echo $run'_toxpoly.bb'`
outfile=${run}.toxpoly.ss
runMB $run $outfile $bb

bb=`echo $run'_toxai.bb'`
outfile=${run}.toxai.ss
runMB $run $outfile $bb

bb=`echo $run'_toxsa.bb'`
outfile=${run}.toxsa.ss
runMB $run $outfile $bb

bb=`echo $run'_toxsi.bb'`
outfile=${run}.toxsi.ss
runMB $run $outfile $bb
