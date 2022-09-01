#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 5           	# Request 20 cores
#$ -l h_rt=240:00:00    # Request max hours runtime
#$ -l h_vmem=5G
#$ -l highmem
#$ -j y			# Combine job outputs into 1 log file

LINEAGE=$(awk '{print $1}' lineages | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate raxml-ng

for MARKER in LSU
do

	if [[ -f "alignments/${LINEAGE}/${MARKER}.fasta_genetree" ]]
	then

		#Convert file for RAxML-NG

		raxml-ng --parse \
		         --msa alignments/${LINEAGE}/${MARKER}_alntrimmed.fasta_genetree \
        		 --model GTR+G \
		         --prefix raxmlng/${LINEAGE}/${LINEAGE}_${MARKER}

		#Run ML tree search and bootstrapping for 1000 iterations

		raxml-ng --all \
        		 --msa alignments/${LINEAGE}/${MARKER}_alntrimmed.fasta_genetree \
	        	 --model GTR+G \
		         --prefix raxmlng/${LINEAGE}/${LINEAGE}_${MARKER} \
        		 --seed 2 \
	        	 --threads ${NSLOTS} \
	        	 --bs-trees autoMRE{1000}

	fi

done
