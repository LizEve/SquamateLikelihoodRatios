#!/bin/bash

gblocks_path="/home/gmount/Gblocks_0.91b/Gblocks"

for filename in *.aln;
do
x=`grep -c "^>" $filename`
b1=$((1+(($x*50)/100)))
b2=$((($x*85)/100))
if (($b2<$b1)); then
b2=$b1
fi
$gblocks_path $filename -t=d -b1=$b1 -b2=$b2 -b3=8 -b4=10 -b5=h -p=n
done
