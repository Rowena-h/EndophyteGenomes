#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 4 cors
#$ -l h_rt=0:10:0 	# Request 30 min runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load anaconda3
conda activate quast

quast.py 	../denovo_assembly/abyss/${STRAIN}/k56/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k64/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k72/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k80/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k88/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k96/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k104/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k112/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k120/${STRAIN}-contigs.fa \
		../denovo_assembly/abyss/${STRAIN}/k128/${STRAIN}-contigs.fa \
		-o ${STRAIN}/${STRAIN}_abyss_quast_results \
		-t ${NSLOTS} --fungus
