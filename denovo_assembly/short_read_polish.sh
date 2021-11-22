#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 4 cores
#$ -l h_rt=48:0:0 	# Request 48 hours runtime
#$ -l h_vmem=20G   	# Request 5GB RAM
#$ -l highmem
#$ -j y
#$ -m bea

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains_shortread)

module load bwa
module load samtools
module load seqtk
module load anaconda3
conda activate pilon

for ASSEMBLER in abyss megahit spades
do
	mkdir ${ASSEMBLER}/${STRAIN}/bwa

	#Index assembly	
	bwa index ${ASSEMBLER}/${STRAIN}/${STRAIN}-contigs.fa
	#Align reads to assembly and sort by readname
	bwa mem ${ASSEMBLER}/${STRAIN}/${STRAIN}-contigs.fa ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz -t ${NSLOTS} | samtools sort -@ ${NSLOTS} -n -o ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapped_sorted.bam -
	#Coordinate sort and mark duplicates
	samtools fixmate -m -@ ${NSLOTS} ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapped_sorted.bam - | samtools sort -@ ${NSLOTS} - | samtools markdup -@ ${NSLOTS} - ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapped_sorted_dups.bam
	#Calculate statistics
	samtools flagstat ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapped_sorted_dups.bam > ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapstats
	
	#Coordinate sort for polishing
	samtools sort -@ ${NSLOTS} ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapped_sorted.bam -o ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapped_coordinatesorted.bam
	#Index
	samtools index ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapped_coordinatesorted.bam

	#Polish with pilon
	pilon --genome ${ASSEMBLER}/${STRAIN}/${STRAIN}-contigs.fa --frags ${ASSEMBLER}/${STRAIN}/bwa/${STRAIN}_${ASSEMBLER}_mapped_coordinatesorted.bam --output ${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished --changes --fix all

	#Remove contigs <200bp (to be NCBI compliant)
	seqtk seq -L 200 ${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished.fasta > ${ASSEMBLER}/${STRAIN}/${STRAIN}_${ASSEMBLER}_polished_filtered.fa
done
