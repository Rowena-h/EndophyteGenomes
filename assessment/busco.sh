#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 4 cors
#$ -l h_rt=48:00:00 	# Request 30 min runtime
#$ -l h_vmem=2G   	# Request 1GB RAM
#$ -m bea
#$ -j y

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load busco

export AUGUSTUS_CONFIG_PATH=/data/SBCS-BuggsLab/RowenaHill/genome_assemblies/augustus_config/config

BUSCO.py -i ../polishing/${STRAIN}/${STRAIN}_abyss_pilon.fasta -c ${NSLOTS} -o ${STRAIN}_abyss_pilon_busco -m genome -l /data/SBCS-BuggsLab/RowenaHill/genome_assemblies/busco_datasets/ascomycota_odb10.2020-09-10

mv run_${STRAIN}_abyss_pilon_busco ${STRAIN}
