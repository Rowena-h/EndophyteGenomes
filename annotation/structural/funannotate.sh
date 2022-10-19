#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=240:00:0 	# Request 240 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -j y

STRAIN=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
ASSEMBLER=$(cat ../../strains_shortread ../../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)

NAME=$(grep ${STRAIN} strain_info | awk -F'\t' '{print $2}')
LOCUS_TAG=$(grep ${STRAIN} strain_info | awk -F'\t' '{print $3}')
TRANSCRIPTS=$(grep ${STRAIN} strain_info | awk -F'\t' '{print $4}')
PROTEINS=$(grep ${STRAIN} strain_info | awk -F'\t' '{print $5}')

#~/funannotate_latest.sif
#funannotate sort	-i ../repeat_masking/${STRAIN}_masked/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_ncbi.fa.masked \
#			-o ../repeat_masking/${STRAIN}_masked/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_ncbi.fa.masked.sorted


#~/funannotate_latest.sif
/data/scratch/btx494/funannotate_v1.8.13.sif funannotate predict	-i ../repeat_masking/${STRAIN}_masked/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_ncbi.fa.masked.sorted \
			--species "${NAME}" \
			--strain IMI${STRAIN} \
			--name ${LOCUS_TAG} \
			--transcript_evidence ${TRANSCRIPTS} \
			--protein_evidence ${PROTEINS} \
			--no-evm-partitions \
			-o ${STRAIN} \
			--cpus ${NSLOTS}
