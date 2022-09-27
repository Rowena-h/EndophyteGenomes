#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 10 		# Request 10 cores
#$ -l h_rt=240:00:00 	# Request 240 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m ea
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)
ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)

mkdir $STRAIN

module load blast+/2.11.0
#module load anaconda3
#conda activate blast

blastn	-query ../denovo_assembly/${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa \
	-db /data/scratch/btx494/nt/nt \
	-task megablast \
	-outfmt '6 qseqid staxids bitscore std' \
	-max_target_seqs 1 \
	-max_hsps 1 \
	-evalue 1e-25 \
	-out ${STRAIN}/${STRAIN}_${ASSEMBLER}_blast.tsv \
	-num_threads ${NSLOTS}
