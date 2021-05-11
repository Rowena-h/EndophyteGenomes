#!/bin/sh
#Script to submit assembly polishing jobs

STRAINS=$(cat ../strains2)
NUM=$(cat ../strains2 | wc -l)

for STRAIN in $STRAINS
do
        mkdir ${STRAIN}
done

qsub -t 1-${NUM} polish.sh
