#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 10 		# Request 10 cores
#$ -l h_rt=240:00:00 	# Request 240 hours runtime
#$ -l highmem
#$ -l h_vmem=60G
#$ -m bea
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)

mkdir $STRAIN

module load anaconda3
conda activate diamond

if grep -Fq ${STRAIN} ../strains_shortread
then

	diamond	blastx \
		--index-chunks 1 \
		-q ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa \
		-d /data/scratch/btx494/uniref90.fasta.dmnd \
		--outfmt 6 \
		-e 1e-25 \
		--out ${STRAIN}/${STRAIN}_spades_diamond.tsv \
		-p ${NSLOTS} \
		--sensitive

elif grep -Fq ${STRAIN} ../strains_hybrid
then

	ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | awk '{print $4}' | sed -n ${SGE_TASK_ID}p)

	diamond blastx \
		--index-chunks 1 \
		-q ../denovo_assembly/${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa \
                -d /data/scratch/btx494/uniref90.fasta.dmnd \
                --outfmt 6 \
                -e 1e-25 \
                --out ${STRAIN}/${STRAIN}_${ASSEMBLER}_diamond.tsv \
                -p ${NSLOTS} \
                --sensitive

fi
