#!/bin/bash
# Move into gene folder
cd $1
WDIR=`pwd`
run=`echo $(basename $WDIR )`

# Currently points to 2.5 version, with jeremy's edits 
function runMB {
    run=$1
    outfile=$2
    bb=$3
    ckp=${outfile%.ss}.ckp
    if [ -f $outfile ] # if outfile exists
    then
        if [[ $(wc -l < "$outfile") -eq 57 ]] # I dont think this happens, either there is an ss file or there isnt. just in case leaving this in. 
        then
            echo $outfile" done"
        else
            sed -i.old 's/mcmcp ngen=2000000/mcmcp append=yes ngen=2000000/g' $bb # restart from ckp file 
            mpirun -np 4 /home/gmount/mb_3.2.5_topoHack_mpi $bb # run mrbayes if output isnt long enough
        fi
    elif [ -f $ckp] # if run has started and ckp file exists
    then 
        sed -i.old 's/mcmcp ngen=2000000/mcmcp append=yes ngen=2000000/g' $bb # restart from ckp file 
        mpirun -np 4 /home/gmount/mb_3.2.5_topoHack_mpi $bb # run mrbayes if output isnt long enough

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

bb=`echo $run'_toxai.bb'`
outfile=${run}.toxai.ss
runMB $run $outfile $bb

bb=`echo $run'_toxsa.bb'`
outfile=${run}.toxsa.ss
runMB $run $outfile $bb

bb=`echo $run'_toxsi.bb'`
outfile=${run}.toxsi.ss
runMB $run $outfile $bb
