#!/bin/sh
#Script to submit assembly quality assessment jobs

STRAINS=$(cat ../strains_shortread ../strains_hybrid)
NUM=$(cat ../strains_shortread ..strains_hybrid | wc -l)

for STRAIN in $STRAINS
do
        mkdir ${STRAIN}
done

qsub -t 1-${NUM} quast.sh
qsub -t 1-${NUM} busco.sh
qsub -t 1-${NUM} diamond.sh
qsub -t 1-${NUM} blast.sh
qsub -t 1-${NUM} read_mapping.sh
