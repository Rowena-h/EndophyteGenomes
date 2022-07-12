#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 4 cores
#$ -l h_rt=72:0:0 	# Request 48 hours runtime
#$ -l highmem
#$ -l h_vmem=60G   	# Request 60GB RAM
#$ -j y
#$ -m bea

STRAIN=$(awk '{print $1}' ../strains_shortread | sed -n ${SGE_TASK_ID}p)

module load spades

spades.py 	-1 ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
		-2 ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
		--careful \
		-o spades/${STRAIN} \
		-t ${NSLOTS}

mv spades/${STRAIN}/contigs.fasta spades/${STRAIN}/${STRAIN}-contigs.fa

module load bwa
module load samtools

#Polishing with short reads

#Index assembly
bwa index spades/${STRAIN}/${STRAIN}-contigs.fa

#Align reads to assembly, coordinate sort and mark duplicates
bwa mem spades/${STRAIN}/${STRAIN}-contigs.fa \
	../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
	../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
	-t ${NSLOTS} | \
	samtools fixmate -m -@ ${NSLOTS} - - | \
	samtools sort -@ ${NSLOTS} - | \
	samtools markdup -@ ${NSLOTS} - spades/${STRAIN}/${STRAIN}_spades_mapped_coordinatesorted.bam

#Calculate statistics
samtools flagstat spades/${STRAIN}/${STRAIN}_spades_mapped_coordinatesorted.bam > spades/${STRAIN}/${STRAIN}_spades_mapstats

#Index for polishing
samtools index spades/${STRAIN}/${STRAIN}_spades_mapped_coordinatesorted.bam

conda activate pilon

#Polish with pilon
pilon --genome spades/${STRAIN}/${STRAIN}-contigs.fa --frags spades/${STRAIN}/${STRAIN}_spades_mapped_coordinatesorted.bam --output spades/${STRAIN}/test --changes --fix all

module load seqtk

#Remove contigs <200bp (to be NCBI compliant)
seqtk seq -L 200 spades/${STRAIN}/${STRAIN}_spades_polished.fasta > spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa
