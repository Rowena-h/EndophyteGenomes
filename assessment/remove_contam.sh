#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 5 cores
#$ -l h_rt=1:00:00     # Request 48 hours runtime
#$ -l h_vmem=1G         # Request 2GB RAM
#$ -m bea
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $1}')
CLASS=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $3}')

module load seqtk

if grep -Fq ${STRAIN} ../strains_shortread
then
	awk 'NR==FNR {A[$1]++;next} $6 in A {print $1}' orders_${CLASS} ${STRAIN}/blobtools/order.${STRAIN}_spades_blobtools.blobDB.table.txt | \
	seqtk subseq ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa - > ${STRAIN}/blobtools/${STRAIN}_spades_polished_filtered_nocontam.fa

elif grep -Fq ${STRAIN} ../strains_hybrid
then
	
	ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')

	awk 'NR==FNR {A[$1]++;next} $6 in A {print $1}' orders_${CLASS} ${STRAIN}/blobtools/order.${STRAIN}_${ASSEMBLER}_blobtools.blobDB.table.txt | \
	seqtk subseq ../denovo_assembly/${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa - > ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa
fi
