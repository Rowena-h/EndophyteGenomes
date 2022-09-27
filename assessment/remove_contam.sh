#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 core
#$ -l h_rt=1:00:00     	# Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $1}')
ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')
CLASS=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $3}')

module load seqtk

#Filter by taxonomy (same class)

awk 'NR==FNR {A[$1]++;next} $6 in A {print $1}' orders_${CLASS} ${STRAIN}/blobtools/order.${STRAIN}_${ASSEMBLER}_blobtools_uniref_nt.blobDB.bestsum.table.txt | \
seqtk subseq ../denovo_assembly/${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa - > ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa

#Filter by coverage (>10)

sed '/^#/d' ${STRAIN}/blobtools/order.${STRAIN}_${ASSEMBLER}_blobtools_uniref_nt.blobDB.bestsum.table.txt | awk '{if($5>10)print $1}' | \
seqtk subseq ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa - > ${STRAIN}/blobtools/${STRAIN}tmp && mv ${STRAIN}/blobtools/${STRAIN}tmp ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa
