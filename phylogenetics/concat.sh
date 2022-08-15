#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 core
#$ -l h_rt=01:00:00 	# Request 1 hour runtime. Max is 240, anything 1 and under has queue priority
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y			# Combine job outputs into 1 log file

LINEAGE=$(awk '{print $1}' lineages | sed -n ${SGE_TASK_ID}p)

#Concatenate gene alignments

module load anaconda3
conda activate AMAS

AMAS.py concat -f fasta -d dna -i alignments/${LINEAGE}/*_alntrimmed.fasta -p alignments/${LINEAGE}/${LINEAGE}_partition.txt -t alignments/${LINEAGE}/${LINEAGE}_concat.phy -u phylip

#Add gene models
sed -i 's/^/GTR+G, /' alignments/${LINEAGE}/${LINEAGE}_partition.txt
