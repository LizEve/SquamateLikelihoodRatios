#!/bin/bash
# Move into gene folder
cd $1
WDIR=`pwd`
run=`echo $(basename $WDIR )`

# Currently points to 2.6, no bug version 
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
            mpirun -np 4 /home/gmount/mb_3.2.5_topoHack_mpi $bb # run mrbayes if output isnt long enough
        fi
    else
    mpirun -np 4 /home/gmount/mb_3.2.5_topoHack_mpi $bb # run mrbayes if output doesnt exists. 
    fi
}


bb=`echo $run'_sclero.bb'`
outfile=${run}.sclero.ss
runMB $run $outfile $bb

bb=`echo $run'_toxpoly.bb'`
outfile=${run}.toxpoly.ss
runMB $run $outfile $bb
