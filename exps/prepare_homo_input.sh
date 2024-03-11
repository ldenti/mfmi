#!/usr/bin/env sh

set -xe

AGC=$1
WD=$2
SEED=23
N=65

mkdir -p $WD
agc listset $AGC | (RANDOM=$SEED; while read line ; do echo "$RANDOM $line" ; done ) | sort | head -n $N | cut -f2 -d' ' > $WD/hlist.txt

# Get sequences of selected haplotypes
mkdir $WD/hfas
while read line
do
    agc getset $AGC $line > $WD/hfas/$line.fa
    samtools faidx $WD/hfas/$line.fa
done < $WD/hlist.txt

# Split by number of haplotypes
for n in 1 2 4 8 16 32 64
do
    mkdir -p $WD/$n/
    head -n $n $WD/hlist.txt > $WD/$n/hlist.txt
    while read idx
    do
	cat $WD/hfas/$idx.fa
    done < $WD/$n/hlist.txt > $WD/$n/reference.fa
    samtools faidx $WD/$n/reference.fa
done

q=$(tail -1 $WD/hlist.txt)
echo $q > $WD/query.txt
ln -sf $WD/hfas/$q.fa $WD/query.fa

rm -rf $WD/hfas
