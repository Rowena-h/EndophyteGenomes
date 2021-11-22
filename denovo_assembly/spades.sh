#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 4 cores
#$ -l h_rt=72:0:0 	# Request 48 hours runtime
#$ -l highmem
#$ -l h_vmem=60G   	# Request 60GB RAM
#$ -j y
#$ -m bea

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains_shortread)

module load spades

spades.py 	-1 ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
		-2 ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
		--careful \
		-o spades/${STRAIN} \
		-t ${NSLOTS}

mv spades/${STRAIN}/contigs.fasta spades/${STRAIN}/${STRAIN}-contigs.fa
