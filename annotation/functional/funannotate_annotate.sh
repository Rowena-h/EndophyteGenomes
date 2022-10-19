#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=240:00:0 	# Request 240 hours runtime
#$ -l h_vmem=60G   	# Request 1GB RAM
#$ -l highmem
#$ -j y
#$ -m ea

STRAIN=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
ASSEMBLER=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)

NAME=$(grep ${STRAIN} ../structural/strain_info | awk -F'\t' '{print $2}')
NAME2=$(grep ${STRAIN} ../structural/strain_info | awk -F'\t' '{print $2}' | sed 's/ /_/g')

~/funannotate_latest.sif funannotate annotate	--gff ../structural/${STRAIN}/predict_results/${NAME2}_IMI${STRAIN}.gff3 \
						--fasta ../structural/${STRAIN}/predict_results/${NAME2}_IMI${STRAIN}.scaffolds.fa \
						--species "${NAME}" \
						--strain IMI${STRAIN} \
						--eggnog eggnog-mapper/${STRAIN}_proteins.emapper.annotations \
						--iprscan interproscan/${STRAIN}_proteins_interproscan.xml \
						--antismash antismash/${STRAIN}/${STRAIN}_proteins.gbk \
						--sbt ../../ncbi_upload/template_${STRAIN}.sbt \
						-o ${STRAIN}_functionalannotation \
						--cpus ${NSLOTS}
