#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=01:00:00
#$ -j y
#$ -t 1

STRAIN=$(cat ../strains_shortread ../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
NAME=$(grep ${STRAIN} ../annotation/funannotate/strain_info | awk -F'\t' '{print $2}')
NAME2=$(grep ${STRAIN} ../annotation/funannotate/strain_info | awk -F'\t' '{print $2}' | sed 's/ /_/g')
ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)
LOCUS_TAG=$(grep ${STRAIN} ../annotation/funannotate/strain_info | awk -F'\t' '{print $3}')

mkdir ${STRAIN}

~/Programmes/table2asn 	-i ../annotation/funannotate/${STRAIN}/predict_results/${NAME2}_IMI${STRAIN}.scaffolds.fa \
			-outdir ${STRAIN} \
			-f ../annotation/funannotate/${STRAIN}//predict_results/${NAME2}_IMI${STRAIN}.gff3 \
			-j "[organism=${NAME}][strain=IMI${STRAIN}]" \
			-t template_${STRAIN}.sbt \
			-locus-tag-prefix ${LOCUS_TAG} \
			-c w \
			-M n \
			-Z \
			-euk

cp ${STRAIN}/${NAME2}_IMI${STRAIN}.scaffolds.sqn aspera_upload/${STRAIN}.sqn
