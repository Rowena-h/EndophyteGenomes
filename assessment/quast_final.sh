#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 4 cors
#$ -l h_rt=0:10:0 	# Request 30 min runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

STRAINS=$(cat ../strains_shortread ../strains_hybrid | awk '{print $1}')
ASSEMBLERS=$(cat ../strains_shortread ../strains_hybrid | awk '{print $4}')

FIRST=$(echo $STRAINS | awk '{print $1}')
FIRST_ASSEMBLER=$(echo $ASSEMBLERS | awk '{print $1}')
LIST=$(echo $FIRST | sed 's#$#/blobtools/'${FIRST}'_'${FIRST_ASSEMBLER}'_polished_filtered_nocontam_ncbi.fa#')

STRAINS2=$(echo $STRAINS | awk '{ $1=""; print}')

for STRAIN in $STRAINS2
do

	ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | grep $STRAIN | awk '{print $4}')

	STRAIN=$(echo $STRAIN | sed 's#$#/blobtools/'${STRAIN}'_'${ASSEMBLER}'_polished_filtered_nocontam_ncbi.fa#')

	LIST=$(echo $LIST $STRAIN)

done

module load anaconda3
conda activate quast

quast.py	${LIST} \
		-l $(echo ${STRAINS} | sed 's/ /,/g') \
                -o final_quast_results \
                -t ${NSLOTS} --fungus
