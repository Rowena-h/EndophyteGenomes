#!/bin/sh

module bedtools
module blast+

for LINEAGE in $(awk '{print $1}' ../lineages}
do
	
	STRAINS=$(grep ${LINEAGE} ../lineages | awk '{print $2}')
	
	for STRAIN in ${STRAINS//,/ }
	do

		ASSEMBLY=$(ls ../../assessment/${STRAIN}/blobtools/${STRAIN}*polished_filtered_nocontam.fa)

		for MARKER in $(cat ../markers)
		do

        		if [[ -f "${LINEAGE}/${MARKER}_example.fasta" ]]
		        then

				#Pull all hits for each marker from the assemblies
				echo "a" | ~/Scripts/GenePull -g ${LINEAGE}/${MARKER}_example.fasta -a ${ASSEMBLY} -o ${LINEAGE}/${STRAIN}_${MARKER}
		
			fi

		done

	done

done
