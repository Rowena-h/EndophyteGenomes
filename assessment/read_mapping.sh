#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 5            # Request 5 cores
#$ -l h_rt=24:00:00     # Request 48 hours runtime
#$ -l h_vmem=2G         # Request 2GB RAM
#$ -m bea
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $1}')

module load bwa
module load samtools
module load minimap2

if grep -Fq ${STRAIN} ../strains_shortread
then
	#Index assembly
	bwa index ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa
	#Align reads to assembly, coordinate sort and mark duplicates
	bwa mem ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa \
		../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
		../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
		-t ${NSLOTS} | \
		samtools fixmate -m -@ ${NSLOTS} - - | \
		samtools sort -@ ${NSLOTS} - | \
		samtools markdup -@ ${NSLOTS} - ${STRAIN}/${STRAIN}_spades_srmapped_coordinatesorted.bam

elif grep -Fq ${STRAIN} ../strains_hybrid
then

	TOOL=$(cat ../strains_shortread ../strains_hybrid | sed -n ${SGE_TASK_ID}p | awk '{print $4}')

	minimap2	-ax map-ont \
			-t ${NSLOTS} \
			-L ../denovo_assembly/${TOOL}/${STRAIN}/${STRAIN}_${TOOL}_polished_filtered.fa \
			../reads/${STRAIN}/${STRAIN}_minion_all_passed.fastq.gz  | \
			samtools fixmate -m -@ ${NSLOTS} - - | \
			samtools sort -@ ${NSLOTS} - | \
			samtools markdup -@ ${NSLOTS} - ${STRAIN}/${STRAIN}_${TOOL}_lrmapped_coordinatesorted.bam

fi
