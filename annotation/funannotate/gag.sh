#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 core
#$ -l h_rt=1:00:0       # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

STRAIN=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
NAME=$(grep ${STRAIN} ../strain_info | awk -F'\t' '{print $2}')
NAME2=$(grep ${STRAIN} ../strain_info | awk -F'\t' '{print $2}' | sed 's/ /_/g')
ASSEMBLER=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)
LOCUS_TAG=$(grep ${STRAIN} ../strain_info | awk -F'\t' '{print $3}')

python2.7 ~/Programmes/genomeannotation-GAG-997e384/gag.py \
-f ../repeat_masking/${STRAIN}_masked/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa.masked.sorted \
-g ${STRAIN}/predict_results/${NAME2}_IMI${STRAIN}.gff3 \
-o ${STRAIN}_gag \
-ris 10 \
--fix_terminal_ns \
--fix_start_stop \

cd ${STRAIN}_gag

rename genome ${STRAIN} *
