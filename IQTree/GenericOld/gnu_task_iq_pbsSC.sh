#!/bin/bash
constraint='scleroglossa.constraint'
# Move into gene folder
cd $1
cp '../'$constraint .
r=$(basename $PWD)
infile=${r}.ntg.nex
prefix=${r}.${constraint%%.*}
echo $prefix
/home/gmount/iqtree -s $infile -pre $constraint -m GTR+G -nt 2 -wsl -g $constraint --runs 10 -b 100 -pers 0.2 -nstop 500 -quiet -redo 