
#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1 		# Request 1 core
#$ -l h_rt=1:00:00 	# Request 1 hour runtime
#$ -l h_vmem=3G   	# Request 3GB RAM
#$ -j y
#$ -m bea

STRAIN=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $1}')

module load samtools
module load anaconda3
conda activate blobtools

mkdir ${STRAIN}/blobtools

if grep -Fq ${STRAIN} ../strains_shortread
then

	samtools index ${STRAIN}/${STRAIN}_spades_srmapped_coordinatesorted.bam

	~/Programmes/blobtools/blobtools create -i ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa \
						-b ${STRAIN}/${STRAIN}_spades_srmapped_coordinatesorted.bam \
						-t ${STRAIN}/${STRAIN}_spades_blast.tsv \
						-o ${STRAIN}/blobtools/${STRAIN}_spades_blobtools
	
	for i in order
	do

		~/Programmes/blobtools/blobtools view  -i ${STRAIN}/blobtools/${STRAIN}_spades_blobtools.blobDB.json \
 	                                               -o ${STRAIN}/blobtools/${i} \
 	                                               --rank $i

		~/Programmes/blobtools/blobtools plot 	-r $i \
							-i ${STRAIN}/blobtools/${STRAIN}_spades_blobtools.blobDB.json \
							-o ${STRAIN}/blobtools/
	done

elif grep -Fq ${STRAIN} ../strains_hybrid
then

	TOOL=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')
	
	samtools index ${STRAIN}/${STRAIN}_${TOOL}_lrmapped_coordinatesorted.bam

        ~/Programmes/blobtools/blobtools create         -i ../denovo_assembly/${TOOL}/${STRAIN}/${STRAIN}_${TOOL}_polished_filtered.fa \
                                                        -b ${STRAIN}/${STRAIN}_${TOOL}_lrmapped_coordinatesorted.bam \
                                                       	-t ${STRAIN}/${STRAIN}_${TOOL}_blast.tsv \
                                                        -o ${STRAIN}/blobtools/${STRAIN}_${TOOL}_blobtools

	for i in order
        do
                ~/Programmes/blobtools/blobtools plot   -r $i \
                                                        -i ${STRAIN}/blobtools/${STRAIN}_${TOOL}_blobtools.blobDB.json \
                                                        -o ${STRAIN}/blobtools/
		
		~/Programmes/blobtools/blobtools view	-i ${STRAIN}/blobtools/${STRAIN}_${TOOL}_blobtools.blobDB.json \
                                                	-o ${STRAIN}/blobtools/${i} \
	                                                --rank $i

        done

fi
