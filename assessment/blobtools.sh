#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1 		# Request 4 cors
#$ -l h_rt=1:00:00 	# Request 30 min runtime
#$ -l h_vmem=3G   	# Request 1GB RAM
#$ -j y
#$ -m bea

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load samtools

samtools index ../polishing/${STRAIN}/${STRAIN}_abyss_mapped_sorted_dups.bam

module load anaconda3
conda activate blobtools

~/Programmes/blobtools/blobtools create 	-i ../denovo_assembly/abyss/${STRAIN}/${STRAIN}-contigs.fa \
						-b ../polishing/${STRAIN}/${STRAIN}_abyss_mapped_sorted_dups.bam \
						-t ${STRAIN}/${STRAIN}_abyss_blast.tsv \
						-o ${STRAIN}_abyss_blobtools

~/Programmes/blobtools/blobtools plot 	-r species \
					-i ${STRAIN}_abyss_blobtools.blobDB.json

~/Programmes/blobtools/blobtools plot   -r family \
                                        -i ${STRAIN}_abyss_blobtools.blobDB.json

~/Programmes/blobtools/blobtools plot   -r phylum \
                                        -i ${STRAIN}_abyss_blobtools.blobDB.json

mv ${STRAIN}_abyss_blobtools* ${STRAIN}
