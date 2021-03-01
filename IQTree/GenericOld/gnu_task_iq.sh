#!/bin/bash
# Move into gene folder
cd $1
r=$(basename $PWD)
infile=${r}.prg.nex

iqtree -s $infile -pre $r -m GTR+G -nt 2 --runs 10 -b 100 -pers 0.2 -nstop 500 -czb -quiet -redo 