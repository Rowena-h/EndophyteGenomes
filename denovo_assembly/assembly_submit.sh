#!/bin/sh
#Script to submit de novo assembly jobs

STRAINS=$(cat ../strains)

mkdir abyss

for STRAIN in $STRAINS
do
	mkdir abyss/${STRAIN}
done

qsub abyss.sh
