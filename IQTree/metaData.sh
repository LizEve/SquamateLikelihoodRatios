#!/bin/bash
outsummary=$1
# Creates header line to specify marginal likelihoods collected
echo -ne "Locus\t"	> $outsummary
echo -ne "Sequences\t"	>> $outsummary
echo -ne "Columns\t"	>> $outsummary
echo -ne "Dist_Pat\t"	>> $outsummary
echo -ne "Pars_Info\t"	>> $outsummary
echo -ne "Sing_Sites\t"	>> $outsummary
echo -ne "Cons_Sites\t" >> $outsummary
echo -ne "Chi2_Fail\t"	>> $outsummary
echo "Gaps_Ambig" >> $outsummary

#f=Reeder_DNA_ADNP.scleroglossa.constraint.log
#f=Reeder_DNA_AHR.scleroglossa.constraint.log

for f in *log
do

echo $f 

echo -ne ${f%%.*}'\t'	>> $outsummary

seq=`grep -m1 -w "sequences with" $f | awk '{print $3}'`

echo -ne $seq'\t'	>> $outsummary

col=`grep -m1 -w "columns" $f | awk '{print $6}'`

echo -ne $col'\t'	>> $outsummary

dis=`grep -m1 -w "distinct" $f | awk '{print $8}'`

echo -ne $dis'\t'	>> $outsummary

par=`grep -m1 -w "parsimony-informative" $f | awk '{print $1}'`

echo -ne $par'\t'	>> $outsummary

sin=`grep -m1 -w "singleton" $f | awk '{print $3}'`

echo -ne $sin'\t'	>> $outsummary

con=`grep -m1 -w "constant" $f | awk '{print $6}'`

echo -ne $con'\t'	>> $outsummary

chi=`grep -m1 -w "sequences failed" $f | awk '{print $4}'`

echo -ne $chi'\t'	>> $outsummary

gap=`grep -m1 -w "sequences contain" $f | awk '{print $2}'`

if [ -z "$gap" ] 
then
	echo "0"	>> $outsummary
else
	echo $gap	>> $outsummary
fi

done
