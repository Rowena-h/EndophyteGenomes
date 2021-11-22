#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1 		# Request 1 core
#$ -l h_rt=1:00:00 	# Request 1 hour runtime
#$ -l h_vmem=3G   	# Request 3GB RAM
#$ -j y
#$ -m bea

STRAIN=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p)

module load samtools
module load anaconda3
conda activate blobtools

mkdir ${STRAIN}/blobtools

if grep -Fq ${STRAIN} ../strains_shortread
then

	samtools index ../spades/${STRAIN}/bwa/${STRAIN}_spades_mapped_sorted_dups.bam

	~/Programmes/blobtools/blobtools create -i ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa \
						-b ../denovo_assembly/spades/${STRAIN}/bwa/${STRAIN}_spades_mapped_sorted_dups.bam \
						-t ${STRAIN}/${STRAIN}_spades_blast.tsv \
						-o ${STRAIN}/blobtools/${STRAIN}_spades_blobtools

	for i in species family phylum
	do
		~/Programmes/blobtools/blobtools plot 	-r $i \
							-i ${STRAIN}/blobtools/${STRAIN}_spades_blobtools.blobDB.json \
							-o ${STRAIN}/blobtools/${STRAIN}_spades_blobtools
	done

elif grep -Fq ${STRAIN} ../strains_hybrid
then
	
	 samtools index ../raven/${STRAIN}/${STRAIN}_raven_mapped_sorted_dups.bam

        ~/Programmes/blobtools/blobtools create         -i ../denovo_assembly/raven/${STRAIN}/${STRAIN}_raven_polished_filtered.fa
                                                        -b ../denovo_assembly/raven/${STRAIN}/${STRAIN}_raven_mapped_sorted_dups.bam
                                                        -t ${STRAIN}/${STRAIN}_raven_blast.tsv \
                                                        -o ${STRAIN}/blobtools/${STRAIN}_raven_blobtools

	for i in species family phylum
        do
                ~/Programmes/blobtools/blobtools plot   -r $i \
                                                        -i ${STRAIN}/blobtools/${STRAIN}_raven_blobtools.blobDB.json \
                                                        -o ${STRAIN}/blobtools/${STRAIN}_raven_blobtools
        done

fi
