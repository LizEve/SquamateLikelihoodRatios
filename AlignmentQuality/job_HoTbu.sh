#!/bin/bash

dataset='Burbrink'
input='gene.list'

working=`pwd`
storage='/media/ExtraDrive3/HoT/Burbrink_HoT'

# Move files make sure to check paths in python script 
#python3 MoveNexToFasta.py $dataset

# Make list of fasta files
cd $working
cp $storage"/"$input $working
JOBS_PER_NODE=2
TASK=$working"/RunHoT.sh"
FILES=$working"/"$input

# Set up parallel call
# set log file, pass workers, whatever nodefile is, and working directory
PARALLEL="parallel -j $JOBS_PER_NODE --joblog log_$dataset.out --wd $working"

# Now call parallel and pass files and task.
$PARALLEL -a $FILES sh $TASK {} 


