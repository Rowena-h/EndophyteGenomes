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

	BUSCO.py	-i ${STRAIN}/blobtools/${STRAIN}_spades_polished_filtered_nocontam.fa \
			-c ${NSLOTS} \
			-o ${STRAIN}_spades_final_busco \
			-m genome \
			-l ../../busco_datasets/ascomycota_odb10.2020-09-10

elif grep -Fq ${STRAIN} ../strains_hybrid
then
	
	ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')

	BUSCO.py 	-i ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa \
			-c ${NSLOTS} \
			-o ${STRAIN}_${ASSEMBLER}_final_busco \
			-m genome \
			-l ../../busco_datasets/ascomycota_odb10.2020-09-10

fi


mv run_${STRAIN}_*_final_busco/short_summary_${STRAIN}_*_final_busco.txt busco_final
rm -r run_${STRAIN}_*_final_busco
