#!/bin/sh
#Script to submit assembly quality assessment jobs

STRAINS=$(cat ../strains)
NUM=$(cat ../strains | wc -l)

for STRAIN in $STRAINS
do
        mkdir ${STRAIN}
done

#qsub -t 1-${NUM} blast.sh
qsub -t 1-${NUM} quast.sh
#qsub -t 1-${NUM} busco.sh

#blobtools.sh to be submitted after blast.sh has finished
