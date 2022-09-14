#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=01:00:00
#$ -j y
#$ -t 1

STRAIN=$(cat ../strains_shortread ../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
NAME=$(grep ${STRAIN} ../annotation/strain_info | awk -F'\t' '{print $2}')
NAME2=$(grep ${STRAIN} ../annotation/strain_info | awk -F'\t' '{print $2}' | sed 's/ /_/g')
ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)
LOCUS_TAG=$(grep ${STRAIN} ../annotation/strain_info | awk -F'\t' '{print $3}')

mkdir ${STRAIN}

~/Programmes/table2asn 	-i ../annotation/funannotate/${STRAIN}_gag/${STRAIN}.fasta \
			-outdir ${STRAIN} \
			-f ../annotation/funannotate/${STRAIN}_gag/${STRAIN}.gff \
			-j "[organism=${NAME}][strain=IMI${STRAIN}]" \
			-t template_${STRAIN}.sbt \
			-locus-tag-prefix ${LOCUS_TAG} \
			-M n \
			-Z \
			-euk

cp ${STRAIN}/${STRAIN}.sqn aspera_upload
