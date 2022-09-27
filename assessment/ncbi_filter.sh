#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 5 cores
#$ -l h_rt=1:00:00     # Request 48 hours runtime
#$ -l h_vmem=1G         # Request 2GB RAM
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $1}')
ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')

cp ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_ncbi.fa

module load bedtools

#Remove sequences flagged as contaminants by NCBI

if [[ -f "${STRAIN}_ncbi_remove.txt" ]]
then

	awk 'NR==FNR{a[$0];next} /^>/{f=($0 in a ? 1 : 0)} !f' ${STRAIN}_ncbi_remove.txt ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_ncbi.fa > ${STRAIN}tmp && mv ${STRAIN}tmp ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_ncbi.fa

fi

#Trim regions flagged as adapter/mitochondrial contaminants by NCBI

if [[ -f "${STRAIN}_ncbi_trim.bed" ]]
then

	bedtools maskfasta -fi ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_ncbi.fa -bed ${STRAIN}_ncbi_trim.bed -fo ${STRAIN}tmp.fa -mc X

	sed 's/X//g' ${STRAIN}tmp.fa > ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_ncbi.fa

	rm ${STRAIN}tmp.fa

fi
