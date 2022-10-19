#!/bin/sh
#$ -cwd
#$ -pe parallel 48
#$ -l infiniband=sdv-i
#$ -l h_rt=240:0:0
#$ -j y
#$ -m ea

STRAIN=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
NAME=$(grep ${STRAIN} ../structural/strain_info | awk -F'\t' '{print $2}' | sed 's/ /_/g')

module load java/11.0.2

/data/scratch/btx494/interproscan-5.59-91.0/interproscan.sh \
	-i /data/SBCS-BuggsLab/RowenaHill/genome_assemblies/assembly_pipeline/annotation/structural/${STRAIN}/predict_results/${NAME}_IMI${STRAIN}.proteins.fa \
	-o interproscan/${STRAIN}_proteins_interproscan.xml \
	-appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,Phobius,PIRSF,PRINTS,SignalP_EUK,SFLD,SMART,SUPERFAMILY,TIGRFAM \
	-f xml \
	-goterms \
	-pa \
	-cpu ${NSLOTS}
