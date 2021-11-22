#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 4 cors
#$ -l h_rt=0:10:0 	# Request 30 min runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

STRAIN=$(cat ../strains_shortread ../strains_hybrid | awk '{print $1}' | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate quast

if grep -Fq ${STRAIN} ../strains_shortread
then

	quast.py        ../denovo_assembly/abyss/${STRAIN}/${STRAIN}_abyss_polished_filtered.fa \
	                ../denovo_assembly/megahit/${STRAIN}/${STRAIN}_megahit_polished_filtered.fa \
        	        ../denovo_assembly/spades/${STRAIN}/${STRAIN}_spades_polished_filtered.fa \
               		-l "ABySS v2.0.2, MEGAHIT v1.2.9, SPAdes v3.11.1" \
	                -o ${STRAIN}/${STRAIN}_shortread_quast_results \
	                -t ${NSLOTS} --fungus

elif grep -Fq ${STRAIN} ../strains_hybrid
then

	quast.py	../denovo_assembly/flye/${STRAIN}/${STRAIN}_flye_polished_filtered.fa \
			../denovo_assembly/raven/${STRAIN}/${STRAIN}_raven_polished_filtered.fa \
 	               	../denovo_assembly/spades_hybrid/${STRAIN}/${STRAIN}_spades_hybrid_polished_filtered.fa \
                	-l "Flye v2.6, Raven v1.6.1, SPAdes v3.11.1" \
	                -o ${STRAIN}/${STRAIN}_hybrid_quast_results \
        	        -t ${NSLOTS} --fungus

fi
