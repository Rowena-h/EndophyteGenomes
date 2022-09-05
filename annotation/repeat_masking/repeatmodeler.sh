#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=24:00:0 	# Request 24 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y
#$ -m ea

STRAIN=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
ASSEMBLER=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)

mkdir ${STRAIN}

cd ${STRAIN}

module load use.dev
module load repeatmodeler

BuildDatabase -name ${STRAIN} ../../../assessment/${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa

RepeatModeler -database ${STRAIN} -engine ncbi -pa 1 -LTRStruct >& ${STRAIN}_repeats.out
