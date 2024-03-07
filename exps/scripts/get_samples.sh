#!/bin/bash

vcf=$1
n=$(($2+1))
SEED=$n

zgrep "^#CHROM" $vcf | cut -f10- | tr "\t" "\n" | (RANDOM=$SEED; while read line ; do echo "$RANDOM $line" ; done ) | sort | head -n $n | cut -f2 -d' '
