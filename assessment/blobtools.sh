#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1 		# Request 1 core
#$ -l h_rt=1:00:00 	# Request 1 hour runtime
#$ -l h_vmem=60G   	# Request 3GB RAM
#$ -l highmem
#$ -j y
#$ -m ea

STRAIN=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $1}')
ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')

module load samtools
module load anaconda3
conda activate blobtools

mkdir ${STRAIN}/blobtools

if grep -Fq ${STRAIN} ../strains_shortread
then

	samtools index ${STRAIN}/${STRAIN}_${ASSEMBLER}_srmapped_coordinatesorted.bam
	
	~/Programmes/blobtools/blobtools taxify -f ${STRAIN}/${STRAIN}_${ASSEMBLER}_diamond.tsv \
						-m /data/scratch/btx494/uniref/uniref90.fasta.taxlist \
						-s 0 -t 1

	mv ${STRAIN}_${ASSEMBLER}_diamond.tsv.taxified.out ${STRAIN}

	~/Programmes/blobtools/blobtools create -i ../denovo_assembly/${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa \
						-b ${STRAIN}/${STRAIN}_${ASSEMBLER}_srmapped_coordinatesorted.bam \
						-t ${STRAIN}/${STRAIN}_${ASSEMBLER}_diamond.tsv.taxified.out \
						-t ${STRAIN}/${STRAIN}_${ASSEMBLER}_blast.tsv \
						-o ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_blobtools_uniref_nt
	
	for i in order
	do

		~/Programmes/blobtools/blobtools view  -i ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_blobtools_uniref_nt.blobDB.json \
 	                                               -o ${STRAIN}/blobtools/${i} \
 	                                               --rank $i

		~/Programmes/blobtools/blobtools plot 	-r $i \
							-i ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_blobtools_uniref_nt.blobDB.json \
							-o ${STRAIN}/blobtools/
	done

elif grep -Fq ${STRAIN} ../strains_hybrid
then

	ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')
	
	samtools index ${STRAIN}/${STRAIN}_${ASSEMBLER}_lrmapped_coordinatesorted.bam

	~/Programmes/blobtools/blobtools taxify -f ${STRAIN}/${STRAIN}_${ASSEMBLER}_diamond.tsv \
                                                -m /data/scratch/btx494/uniref/uniref90.fasta.taxlist \
                                                -s 0 -t 1

	mv ${STRAIN}_${ASSEMBLER}_diamond.tsv.taxified.out ${STRAIN}

        ~/Programmes/blobtools/blobtools create         -i ../denovo_assembly/${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa \
                                                        -b ${STRAIN}/${STRAIN}_${ASSEMBLER}_lrmapped_coordinatesorted.bam \
                                                       	-t ${STRAIN}/${STRAIN}_${ASSEMBLER}_diamond.tsv.taxified.out \
							-t ${STRAIN}/${STRAIN}_${ASSEMBLER}_blast.tsv \
                                                        -o ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_blobtools_uniref_nt

	for i in order
        do
                ~/Programmes/blobtools/blobtools plot   -r $i \
                                                        -i ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_blobtools_uniref_nt.blobDB.json \
                                                        -o ${STRAIN}/blobtools/
		
		~/Programmes/blobtools/blobtools view	-i ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_blobtools_uniref_nt.blobDB.json \
                                                	-o ${STRAIN}/blobtools/${i} \
	                                                --rank $i

        done

fi
