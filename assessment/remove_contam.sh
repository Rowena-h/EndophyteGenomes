#!/bin/sh

STRAINS=$(cat ../strains_hybrid | awk '{print $1}')

module load seqtk

for STRAIN in $STRAINS
do

	CLASS=$(cat ../strains_shortread ../strains_hybrid | grep $STRAIN | awk '{print $3}')

	if grep -Fq ${STRAIN} ../strains_shortread
	then

		awk 'NR==FNR {A[$1]++;next} $6 in A {print $1}' orders_${CLASS} ${STRAIN}/blobtools/${STRAIN}_spades_blobtools_order.${STRAIN}_spades_blobtools.blobDB.table.txt | seqtk subseq ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa - > ${STRAIN}/blobtools/${STRAIN}_spades_polished_filtered_${CLASS}.fa

	elif grep -Fq ${STRAIN} ../strains_hybrid
	then

		awk 'NR==FNR {A[$1]++;next} $6 in A {print $1}' orders_${CLASS} ${STRAIN}/blobtools/${STRAIN}_raven_blobtools_order.${STRAIN}_raven_blobtools.blobDB.table.txt | seqtk subseq ../denovo_assembly/raven/${STRAIN}/${STRAIN}_raven_polished_filtered.fa - > ${STRAIN}/blobtools/${STRAIN}_raven_polished_filtered_${CLASS}.fa

	fi
	
done



