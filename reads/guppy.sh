#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 12     # Request 12 cores
#$ -l h_rt=72:0:0 # Request 72 hours runtime
#$ -l h_vmem=1G   # Request 1GB RAM
#$ -j y

STRAIN=$(awk '{print $1}' ../strains_hybrid | sed -n ${SGE_TASK_ID}p)

module load singularity

singularity exec /data/containers/nanopore/nanopore-guppy-4.5.3.simg guppy_basecaller \
	-i /data/scratch/btx494/${STRAIN}/fast5/ \
	-s ${STRAIN}/minion/ \
	-c dna_r9.4.1_450bps_fast.cfg \
	--compress_fastq \
	--num_callers 4 --cpu_threads_per_caller 3

cat ${STRAIN}/minion/pass/*.fastq.gz > ${STRAIN}/${STRAIN}_minion_all_passed.fastq.gz

#Calculate number of passed bases
for i in ${STRAIN}/minion/pass/*.fastq.gz; do zcat $i | paste - - - - | cut -f 2 | tr -d '\n'| wc -c; done >> ${STRAIN}_tmp
echo "Number of passed bases: `awk '{ sum += $1 } END { print sum }' ${STRAIN}_tmp`"
rm ${STRAIN}_tmp
