#!/bin/sh
#$ -cwd                # Set the working directory for the job to the current directory
#$ -pe smp 10          # Request 4 cores
#$ -l h_rt=12:00:00    # Request 48 hours runtime
#$ -l h_vmem=5G        # Request 5GB RAM
#$ -j y
#$ -m bea

STRAIN=$(awk '{print $1}' ../strains_shortread | sed -n ${SGE_TASK_ID}p)

module load abyss

abyss-fac abyss/${STRAIN}/k*/${STRAIN}-contigs.fa > abyss/${STRAIN}_kcomp

KMER=$(tail -n +2 abyss/${STRAIN}_kcomp | sort -n -r -k 6 | awk '{print $11}' | sed "s#abyss/${STRAIN}/##" | sed "s#/${STRAIN}-contigs.fa##" | head -n +1)

echo "${KMER} selected" >> abyss/${STRAIN}_kcomp

cp abyss/${STRAIN}/${KMER}/${STRAIN}-contigs.fa abyss/${STRAIN}/

module load bwa
module load samtools

#Polishing with short reads

#Index assembly
bwa index abyss/${STRAIN}/${STRAIN}-contigs.fa

#Align reads to assembly, coordinate sort and mark duplicates
bwa mem abyss/${STRAIN}/${STRAIN}-contigs.fa \
	../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
	../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
	-t ${NSLOTS} | \
	samtools fixmate -m -@ ${NSLOTS} - - | \
	samtools sort -@ ${NSLOTS} - | \
	samtools markdup -@ ${NSLOTS} - abyss/${STRAIN}/${STRAIN}_abyss_mapped_coordinatesorted.bam

#Calculate statistics
samtools flagstat abyss/${STRAIN}/${STRAIN}_abyss_mapped_coordinatesorted.bam > abyss/${STRAIN}/${STRAIN}_abyss_mapstats

#Index for polishing
samtools index abyss/${STRAIN}/${STRAIN}_abyss_mapped_coordinatesorted.bam

module load anaconda3
conda activate pilon

#Polish with pilon
pilon --genome abyss/${STRAIN}/${STRAIN}-contigs.fa --frags abyss/${STRAIN}/${STRAIN}_abyss_mapped_coordinatesorted.bam --output abyss/${STRAIN}/${STRAIN}_abyss_polished --changes --fix all

module load seqtk

#Remove contigs <200bp (to be NCBI compliant)
seqtk seq -L 200 abyss/${STRAIN}/${STRAIN}_abyss_polished.fasta > abyss/${STRAIN}/${STRAIN}_abyss_polished_filtered.fa

