#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 core
#$ -l h_rt=01:00:00 	# Request 1 hour runtime. Max is 240, anything 1 and under has queue priority
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y			# Combine job outputs into 1 log file

LINEAGE=$(awk '{print $1}' lineages | sed -n ${SGE_TASK_ID}p)

for MARKER in $(cat markers)
do
	
	if [[ -f "alignments/${LINEAGE}/${MARKER}.fasta" ]]
	then

		#Create alignment

		module load mafft

		mafft alignments/${LINEAGE}/${MARKER}.fasta > alignments/${LINEAGE}/${MARKER}_aln.fasta

		#Trim alignment

		module load anaconda3
		conda activate trimal

		trimal -in alignments/${LINEAGE}/${MARKER}_aln.fasta -fasta -gappyout > alignments/${LINEAGE}/${MARKER}_alntrimmed.fasta

	fi

done
