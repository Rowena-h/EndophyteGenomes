#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 10 		# Request 10 cores
#$ -l h_rt=240:00:00 	# Request 240 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)

mkdir $STRAIN

module load blast+/2.11.0
module load anaconda3
conda activate blast

if grep -Fq ${STRAIN} ../strains_shortread
then

	blastn 	-query ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa \
		-db /data/scratch/btx494/nt \
		-task megablast \
		-outfmt '6 qseqid staxids bitscore std' \
		-max_target_seqs 1 \
		-max_hsps 1 \
		-evalue 1e-25 \
		-out ${STRAIN}/${STRAIN}_spades_blast.tsv \
		-num_threads ${NSLOTS}

elif grep -Fq ${STRAIN} ../strains_hybrid
then

	TOOL=$(cat ../strains_shortread ../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)

	blastn  -query ../denovo_assembly/${TOOL}/${STRAIN}/${STRAIN}_${TOOL}_polished_filtered.fa \
        	-db /data/scratch/btx494/nt \
		-task megablast \
	        -outfmt '6 qseqid staxids bitscore std' \
	        -max_target_seqs 1 \
	        -max_hsps 1 \
	        -evalue 1e-25 \
	        -out ${STRAIN}/${STRAIN}_${TOOL}_blast.tsv \
	        -num_threads ${NSLOTS}

fi
