#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=24:00:00
#$ -l h_vmem=10G
#$ -m bea

STRAIN=$(awk '{print $1}' ../strains_hybrid | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate raven

mkdir raven/${STRAIN}

raven ../reads/${STRAIN}/minion/pass/*.fastq.gz -t ${NSLOTS} > raven/${STRAIN}/${STRAIN}-contigs.fa

conda deactivate

module load bwa
module load samtools

#Polishing with short reads

#Index assembly
bwa index raven/${STRAIN}/${STRAIN}-contigs.fa
#Align reads to assembly and sort by readname
bwa mem raven/${STRAIN}/${STRAIN}-contigs.fa ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz -t ${NSLOTS} | samtools sort -@ ${NSLOTS} -n -o raven/${STRAIN}/${STRAIN}_raven_mapped_sorted.bam -
#Coordinate sort and mark duplicates
samtools fixmate -m -@ ${NSLOTS} raven/${STRAIN}/${STRAIN}_raven_mapped_sorted.bam - | samtools sort -@ ${NSLOTS} - | samtools markdup -@ ${NSLOTS} - raven/${STRAIN}/${STRAIN}_raven_mapped_sorted_dups.bam
#Calculate statistics
samtools flagstat raven/${STRAIN}/${STRAIN}_raven_mapped_sorted_dups.bam > raven/${STRAIN}/${STRAIN}_raven_mapstats

#Coordinate sort for polishing
samtools sort -@ ${NSLOTS} raven/${STRAIN}/${STRAIN}_raven_mapped_sorted.bam -o raven/${STRAIN}/${STRAIN}_raven_mapped_coordinatesorted.bam
#Index
samtools index raven/${STRAIN}/${STRAIN}_raven_mapped_coordinatesorted.bam

conda activate pilon

#Polish with pilon
pilon --genome raven/${STRAIN}/${STRAIN}-contigs.fa --frags raven/${STRAIN}/${STRAIN}_raven_mapped_coordinatesorted.bam --output raven/${STRAIN}/${STRAIN}_raven_polished --changes --fix all

module load seqtk

#Remove contigs <200bp (to be NCBI compliant)
seqtk seq -L 200 raven/${STRAIN}/${STRAIN}_raven_polished.fasta > raven/${STRAIN}/${STRAIN}_raven_polished_filtered.fa
