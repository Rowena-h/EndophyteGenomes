#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 10 		# Request 4 cors
#$ -l h_rt=240:00:00 	# Request 30 min runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -j y

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

#module load blast+/2.6.0
module load anaconda3
conda activate blast

blastn 	-query ../denovo_assembly/abyss/${STRAIN}/${STRAIN}-contigs.fa \
	-db /data/scratch/btx494/nt \
	-outfmt '6 qseqid staxids bitscore std' \
	-max_target_seqs 1 \
	-max_hsps 1 \
	-evalue 1e-25 \
	-out ${STRAIN}/${STRAIN}_abyss_blast.tsv \
	-num_threads ${NSLOTS}
