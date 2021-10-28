#!/bin/sh

SR_STRAINS=$(cat ../strains_shortread)

module load abyss

for STRAIN in $(SR_STRAINS)
do
	KMER=

	abyss-fac abyss/${STRAIN}/k*/${STRAIN}-contigs.fa > abyss/${STRAIN}_kcomp

	KMER=$(tail -n +2 abyss/${STRAIN}_kcomp | sort -r -k 6 | awk '{print $11}' | sed "s#abyss/${STRAIN}/##" | sed "s#/${STRAIN}-contigs.fa##" | head -n +1)

	echo "${KMER} selected" >> abyss/${STRAIN}_kcomp

	cp abyss/${STRAIN}/${KMER}/${STRAIN}-contigs.fa abyss/${STRAIN}/
done
