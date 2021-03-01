#!/bin/bash
# Move into gene folder
cd $1
WDIR=`pwd`
run=`echo $(basename $WDIR )`
bb=`echo $run'.bb'`

/home/gmount/MrBayes/src/mb $bb