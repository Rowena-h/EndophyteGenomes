#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=12:00:0 	# Request 240 hours runtime
#$ -l h_vmem=5G   	# Request 1GB RAM
#$ -j y
#$ -m ea

STRAIN=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
NAME=$(grep ${STRAIN} ../structural/strain_info | awk -F'\t' '{print $2}' | sed 's/ /_/g')

module load anaconda3
conda activate antismash

antismash ../structural/${STRAIN}/predict_results/${NAME}_IMI${STRAIN}.scaffolds.fa \
	--output-dir antismash/${STRAIN} \
	--taxon fungi \
	--cassis \
	--cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go \
	--genefinding-gff3 ../structural/${STRAIN}/predict_results/${NAME}_IMI${STRAIN}.gff3 \
	--output-basename ${STRAIN}_proteins \
	-v \
	-c ${NSLOTS}
