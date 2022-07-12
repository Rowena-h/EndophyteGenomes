#!/bin/sh
#$ -cwd
#$ -pe smp 10
#$ -l h_rt=48:0:0
#$ -j y
#$ -m bea

STRAIN=$(awk '{print $1}' ../strains_shortread | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate megahit

megahit -1 ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
	-2 ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
	-t ${NSLOTS} \
	-o megahit/${STRAIN} \
	--k-list 51,59,67,75,83,91,99,107,115,123,131

mv megahit/${STRAIN}/final.contigs.fa megahit/${STRAIN}/${STRAIN}-contigs.fa

module load bwa
module load samtools

#Polishing with short reads

#Index assembly
bwa index megahit/${STRAIN}/${STRAIN}-contigs.fa

#Align reads to assembly, coordinate sort and mark duplicates
bwa mem megahit/${STRAIN}/${STRAIN}-contigs.fa \
	../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz \
	../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz \
	-t ${NSLOTS} | \
	samtools fixmate -m -@ ${NSLOTS} - - | \
	samtools sort -@ ${NSLOTS} - | \
	samtools markdup -@ ${NSLOTS} - megahit/${STRAIN}/${STRAIN}_megahit_mapped_coordinatesorted.bam

#Calculate statistics
samtools flagstat megahit/${STRAIN}/${STRAIN}_megahit_mapped_coordinatesorted.bam > megahit/${STRAIN}/${STRAIN}_megahit_mapstats

#Index for polishing
samtools index megahit/${STRAIN}/${STRAIN}_megahit_mapped_coordinatesorted.bam

conda activate pilon

#Polish with pilon
pilon --genome megahit/${STRAIN}/${STRAIN}-contigs.fa --frags megahit/${STRAIN}/${STRAIN}_megahit_mapped_coordinatesorted.bam --output megahit/${STRAIN}/test --changes --fix all

module load seqtk

#Remove contigs <200bp (to be NCBI compliant)
seqtk seq -L 200 megahit/${STRAIN}/${STRAIN}_megahit_polished.fasta > megahit/${STRAIN}/${STRAIN}_megahit_polished_filtered.fa
