#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 5            # Request 5 cores
#$ -l h_rt=12:00:00     # Request 48 hours runtime
#$ -l h_vmem=5G         # Request 2GB RAM
#$ -l highmem
#$ -m bea
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $1}')

module load bwa
module load samtools
module load minimap2

if grep -Fq ${STRAIN} ../strains_shortread
then
	#Index assembly
	bwa index ${STRAIN}/blobtools/${STRAIN}_spades_polished_filtered_nocontam.fa
	#Align reads to assembly, coordinate sort, mark duplicates and produce read mapping stats
	bwa mem ${STRAIN}/blobtools/${STRAIN}_spades_polished_filtered_nocontam.fa \
		../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
		../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
		-t ${NSLOTS} | \
		samtools fixmate -m -@ ${NSLOTS} - - | \
		samtools sort -@ ${NSLOTS} - | \
		samtools markdup -@ ${NSLOTS} - - | \
		samtools stats - > ${STRAIN}/blobtools/${STRAIN}_spades_polished_filtered_nocontam_mapstats

	rm ${STRAIN}/blobtools/${STRAIN}_spades_polished_filtered_nocontam.fa.*

elif grep -Fq ${STRAIN} ../strains_hybrid
then

	ASSEMBLER=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')

	#Index assembly
        bwa index ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa
        #Align reads to assembly, coordinate sort, mark duplicates and produce read mapping stats
        bwa mem ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa \
                ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
                ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
                -t ${NSLOTS} | \
                samtools fixmate -m -@ ${NSLOTS} - - | \
                samtools sort -@ ${NSLOTS} - | \
                samtools markdup -@ ${NSLOTS} - - | \
                samtools stats - > ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_srmapstats

	minimap2	-ax map-ont \
			-t ${NSLOTS} \
        	        -L ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam.fa \
               		../reads/${STRAIN}/${STRAIN}_minion_all_passed.fastq.gz  | \
	                samtools fixmate -m -@ ${NSLOTS} - - | \
			samtools sort -@ ${NSLOTS} - | \
			samtools markdup -@ ${NSLOTS} - - | \
			samtools stats - > ${STRAIN}/blobtools/${STRAIN}_${ASSEMBLER}_polished_filtered_nocontam_lrmapstats

fi
