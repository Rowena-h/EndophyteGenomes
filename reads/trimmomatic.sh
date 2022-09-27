#!/bin/sh
#$ -cwd           
#$ -pe smp 4      
#$ -l h_rt=02:00:00 
#$ -l h_vmem=5G   
#$ -j y

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains_shortread | awk '{print $1}')

module load trimmomatic

java -jar /share/apps/centos7/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads ${NSLOTS} -trimlog ${STRAIN}/${STRAIN}-trimmomatic.log ${STRAIN}/${STRAIN}_1.fastq.gz ${STRAIN}/${STRAIN}_2.fastq.gz ${STRAIN}/${STRAIN}_1_trimmedpaired.fastq.gz ${STRAIN}/${STRAIN}_1_trimmedunpaired.fastq.gz ${STRAIN}/${STRAIN}_2_trimmedpaired.fastq.gz ${STRAIN}/${STRAIN}_2_trimmedunpaired.fastq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
