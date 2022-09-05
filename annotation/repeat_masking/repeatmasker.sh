#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=24:00:0 	# Request 24 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m ea
#$ -j y

STRAIN=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
ASSEMBLER=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate repeatmasker

mkdir ${STRAIN}_masked

RepeatMasker 	-e ncbi -lib ${STRAIN}/RM*/consensi.fa \
		-pa ${NSLOTS} -xsmall \
		-dir ${STRAIN}_masked ../../assessment/${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa
