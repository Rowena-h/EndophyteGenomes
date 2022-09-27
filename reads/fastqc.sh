#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 4 cores
#$ -l h_rt=1:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 1GB RAM
#$ -j y

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains_shortread | awk '{print $1}')

module load fastqc

fastqc -t ${NSLOTS} ${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz ${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz
