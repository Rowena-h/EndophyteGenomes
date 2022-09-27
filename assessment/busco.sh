#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=48:00:00 	# Request 48 hours runtime
#$ -l h_vmem=2G   	# Request 2GB RAM
#$ -m bea
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)

module load busco

export AUGUSTUS_CONFIG_PATH=/data/SBCS-BuggsLab/RowenaHill/genome_assemblies/augustus_config/config

if grep -Fq ${STRAIN} ../strains_shortread
then

	for ASSEMBLER in abyss megahit spades
	do
		BUSCO.py 	-i ../denovo_assembly/${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa \
				-c ${NSLOTS} \
				-o ${STRAIN}_${ASSEMBLER}_busco \
				-m genome \
				-l ../../busco_datasets/ascomycota_odb10.2020-09-10

		mv run_${STRAIN}_${ASSEMBLER}* ${STRAIN}

	done

elif grep -Fq ${STRAIN} ../strains_hybrid
then

	for ASSEMBLER in flye raven spades_hybrid
	do
		BUSCO.py 	-i ../denovo_assembly/${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa \
				-c ${NSLOTS} \
				-o ${STRAIN}_${ASSEMBLER}_busco \
				-m genome \
				-l ../../busco_datasets/ascomycota_odb10.2020-09-10

		mv run_${STRAIN}_${ASSEMBLER}* ${STRAIN}

	done

fi
