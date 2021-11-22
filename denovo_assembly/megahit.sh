#!/bin/sh
#$ -cwd
#$ -pe smp 10
#$ -l h_rt=48:0:0
#$ -j y
#$ -m bea

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains_shortread)

module load anaconda3
conda activate megahit

megahit -1 ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
	-2 ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
	-t ${NSLOTS} \
	-o megahit/${STRAIN} \
	--k-list 51,59,67,75,83,91,99,107,115,123,131

mv megahit/${STRAIN}/final.contigs.fa megahit/${STRAIN}/${STRAIN}-contigs.fa
