#!/bin/bash
# $1 is the input from the list that is passed to gnu parallel in the pbs script 
# Move into gene folder
cd $1
r=$(basename $PWD)
infile=${r}.nex

function runIQ {
    infile=$1
    constraint=$2
    #prefix=${constraint%.*}
    outfile1=${constraint}.treefile
    outfile2=${constraint}.iqtree
    # check if both output files exists, they should only be written at end of run i think
    if [[ -f $outfile1 && -f $outfile2 ]]
    then
        echo $outfile1" and "$outfile2" done"
    else
        /home/gmount/iqtree -s $infile -pre $constraint -m GTR+G -nt 1 -wsl -g $constraint --runs 5 -pers 0.2 -nstop 500 -quiet -cptime 20 
    fi
}


constraint=${r}.scleroglossa.constraint
runIQ $infile $constraint
echo $constraint

constraint=${r}.toxicoferaP.constraint
runIQ $infile $constraint
echo $constraint

constraint=${r}.toxicoferaAI.constraint
runIQ $infile $constraint
echo $constraint

constraint=${r}.toxicoferaSA.constraint
runIQ $infile $constraint
echo $constraint

constraint=${r}.toxicoferaSI.constraint
runIQ $infile $constraint
echo $constraint