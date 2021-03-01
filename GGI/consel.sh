#!/bin/bash

FILE=$1
# File must have a .sitelh extention 

CDIR='/home/gmount/consel/bin'
BASE=`basename ${FILE}`
LOCUS="${BASE%.*}"

# create a directory for the output if it doesn't already exist
CONSEL="consel_output"
mkdir -p $CONSEL

# run CONSEL on the file
$CDIR/seqmt --puzzle $LOCUS
$CDIR/makermt $LOCUS
$CDIR/consel $LOCUS
rm -f $LOCUS.rmt
$CDIR/catpv $LOCUS > $LOCUS.out

# move all the output to the CONSEL output directory
mv $LOCUS.* $CONSEL



