#!/bin/sh
#$ -cwd
#$ -pe parallel 48
#$ -l infiniband=sdv-i
#$ -l h_rt=24:0:0
#$ -j y
#$ -t 56-128:8	#kmer sizes to check

STRAINS=$(cat ../strains_shortread)

module load abyss

for STRAIN in $STRAINS
do
	mkdir abyss/${STRAIN}/k${SGE_TASK_ID}
	abyss-pe -C abyss/${STRAIN}/k${SGE_TASK_ID} name=${STRAIN} k=${SGE_TASK_ID} in="/data/SBCS-BuggsLab/RowenaHill/genome_assemblies/assembly_pipeline/reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz /data/SBCS-BuggsLab/RowenaHill/genome_assemblies/assembly_pipeline/reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz" np=${NSLOTS}
done
