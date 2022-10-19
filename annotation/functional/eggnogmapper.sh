#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=12:00:0 	# Request 240 hours runtime
#$ -l h_vmem=5G   	# Request 1GB RAM
#$ -j y
#$ -m ea

STRAIN=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
ASSEMBLER=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)

NAME=$(grep ${STRAIN} ../structural/strain_info | awk -F'\t' '{print $2}')
NAME2=$(grep ${STRAIN} ../structural/strain_info | awk -F'\t' '{print $2}' | sed 's/ /_/g')

module load anaconda3
conda activate eggnog-mapper

emapper.py 	-i ../structural/${STRAIN}/predict_results/${NAME2}_IMI${STRAIN}.proteins.fa \
		-m diamond \
		-o ${STRAIN}_proteins \
		--output_dir eggnog-mapper \
		--cpu ${NSLOTS}
