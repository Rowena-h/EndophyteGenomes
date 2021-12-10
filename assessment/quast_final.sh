#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 4 cors
#$ -l h_rt=0:10:0 	# Request 30 min runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

SR_STRAINS=$(awk '{print $1}' ../strains_shortread)

FIRST=$(echo $SR_STRAINS | awk '{print $1}')
SR_LIST=$(echo $FIRST | sed 's#$#/blobtools/'${FIRST}'_spades_polished_filtered_nocontam.fa#')

SR_STRAINS2=$(echo $SR_STRAINS | awk '{ $1=""; print}')

for STRAIN in $SR_STRAINS2
do

	STRAIN=$(echo $STRAIN | sed 's#$#/blobtools/'${STRAIN}'_spades_polished_filtered_nocontam.fa#')

	SR_LIST=$(echo $SR_LIST $STRAIN)

done

module load anaconda3
conda activate quast

quast.py	${SR_LIST} \
		-l $(echo ${SR_STRAINS} | sed 's/ /,/g') \
                -o shortread_final_quast_results \
                -t ${NSLOTS} --fungus
