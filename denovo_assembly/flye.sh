#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=72:0:0
#$ -l h_vmem=10G
#$ -m bea

STRAIN=$(awk '{print $1}' ../strains_hybrid | sed -n ${SGE_TASK_ID}p)
GENOMESIZE=$(awk '{print $2}' ../strains_hybrid | sed -n ${SGE_TASK_ID}p)

module load flye/2.6

#Long read assembly

flye --nano-raw ../reads/${STRAIN}/minion/pass/*.fastq.gz \
     -g ${GENOMESIZE}m \
     -o flye/${STRAIN} \
     -t ${NSLOTS}

mv flye/${STRAIN}/assembly.fasta flye/${STRAIN}/${STRAIN}-contigs.fa

module unload flye/2.6

module load bwa
module load samtools

#Polishing with short reads

#Index assembly
bwa index flye/${STRAIN}/${STRAIN}-contigs.fa
#Align reads to assembly and sort by readname
bwa mem flye/${STRAIN}/${STRAIN}-contigs.fa ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz -t ${NSLOTS} | samtools sort -@ ${NSLOTS} -n -o flye/${STRAIN}/${STRAIN}_flye_mapped_sorted.bam -
#Coordinate sort and mark duplicates
samtools fixmate -m -@ ${NSLOTS} flye/${STRAIN}/${STRAIN}_flye_mapped_sorted.bam - | samtools sort -@ ${NSLOTS} - | samtools markdup -@ ${NSLOTS} - flye/${STRAIN}/${STRAIN}_flye_mapped_sorted_dups.bam
#Calculate statistics
samtools flagstat flye/${STRAIN}/${STRAIN}_flye_mapped_sorted_dups.bam > flye/${STRAIN}/${STRAIN}_flye_mapstats
	
#Coordinate sort for polishing
samtools sort -@ ${NSLOTS} flye/${STRAIN}/${STRAIN}_flye_mapped_sorted.bam -o flye/${STRAIN}/${STRAIN}_flye_mapped_coordinatesorted.bam
#Index
samtools index flye/${STRAIN}/${STRAIN}_flye_mapped_coordinatesorted.bam

module load anaconda3
conda activate pilon

#Polish with pilon
pilon --genome flye/${STRAIN}/${STRAIN}-contigs.fa --frags flye/${STRAIN}/${STRAIN}_flye_mapped_coordinatesorted.bam --output flye/${STRAIN}/${STRAIN}_flye_polished --changes --fix all

module load seqtk

#Remove contigs <200bp (to be NCBI compliant)
seqtk seq -L 200 flye/${STRAIN}/${STRAIN}_flye_polished.fasta > flye/${STRAIN}/${STRAIN}_flye_polished_filtered.fa
