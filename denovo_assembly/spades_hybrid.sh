#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 4 cores
#$ -l h_rt=72:0:0 	# Request 48 hours runtime
#$ -l highmem
#$ -l h_vmem=60G   	# Request 60GB RAM
#$ -j y
#$ -m bea

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains_hybrid | awk '{print $1}')

module load spades

spades.py 	-1 ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
		-2 ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
		--nanopore ../reads/${STRAIN}/${STRAIN}_minion_all_passed.fastq.gz \
		--careful \
		-o spades_hybrid/${STRAIN} \
		-t ${NSLOTS}

mv spades_hybrid/${STRAIN}/contigs.fasta spades_hybrid/${STRAIN}/${STRAIN}-contigs.fa

module load bwa
module load samtools

#Index assembly
bwa index spades_hybrid/${STRAIN}/${STRAIN}-contigs.fa
#Align reads to assembly and sort by readname
bwa mem spades_hybrid/${STRAIN}/${STRAIN}-contigs.fa ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz -t ${NSLOTS} | samtools sort -@ ${NSLOTS} -n -o spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapped_sorted.bam -
#Coordinate sort and mark duplicates
samtools fixmate -m -@ ${NSLOTS} spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapped_sorted.bam - | samtools sort -@ ${NSLOTS} - | samtools markdup -@ ${NSLOTS} - spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapped_sorted_dups.bam
#Calculate statistics
samtools flagstat spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapped_sorted_dups.bam > spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapstats

#Coordinate sort for polishing
samtools sort -@ ${NSLOTS} spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapped_sorted.bam -o spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapped_coordinatesorted.bam
#Index
samtools index spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapped_coordinatesorted.bam

module load anaconda3
conda activate pilon

#Polish with pilon
pilon --genome spades_hybrid/${STRAIN}/${STRAIN}-contigs.fa --frags spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_mapped_coordinatesorted.bam --output spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_polished --changes --fix all --threads ${NSLOTS}

module load seqtk

#Remove contigs <200bp (to be NCBI compliant)
seqtk seq -L 200 spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_polished.fasta > spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_polished_filtered.fa
