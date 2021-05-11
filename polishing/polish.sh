#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 1 core
#$ -l h_rt=48:0:0 	# Request 24 hour runtime
#$ -l h_vmem=5G   	# Request 1GB RAM
#$ -j y
#$ -m bea

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains2)

module load bwa
module load samtools

#Index assembly	
bwa index ../denovo_assembly/abyss/${STRAIN}/${STRAIN}-contigs.fa
#Align reads to assembly and sort by readname
bwa mem ../denovo_assembly/abyss/${STRAIN}/${STRAIN}-contigs.fa ../reads/${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz ../reads/${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz -t ${NSLOTS} | samtools sort -@ ${NSLOTS} -n -o ${STRAIN}_abyss_mapped_sorted.bam -
#Coordinate sort and mark duplicates
samtools fixmate -m -@ ${NSLOTS} ${STRAIN}_abyss_mapped_sorted.bam - | samtools sort -@ ${NSLOTS} - | samtools markdup -@ ${NSLOTS} - ${STRAIN}_abyss_mapped_sorted_dups.bam
#Calculate statistics
samtools flagstat ${STRAIN}_abyss_mapped_sorted_dups.bam > ${STRAIN}_abyss_mapstats
	
#Coordinate sort for polishing
samtools sort -@ ${NSLOTS} ${STRAIN}_abyss_mapped_sorted.bam -o ${STRAIN}_abyss_mapped_coordinatesorted.bam
#Index
samtools index ${STRAIN}_abyss_mapped_coordinatesorted.bam

#Polish with pilon
java -jar /data/home/btx494/Programmes/pilon/pilon-1.23.jar --genome ../denovo_assembly/abyss/${STRAIN}/${STRAIN}-contigs.fa --frags ${STRAIN}_abyss_mapped_coordinatesorted.bam --output ${STRAIN}_abyss_pilon --changes --fix all --threads ${NSLOTS}

mv ${STRAIN}_abyss* ${STRAIN}
