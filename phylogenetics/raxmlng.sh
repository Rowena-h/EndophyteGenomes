#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 20           # Request 20 cores
#$ -l h_rt=240:00:00    # Request max hours runtime
#$ -l h_vmem=5G
#$ -l highmem
#$ -j y			# Combine job outputs into 1 log file

LINEAGE=$(awk '{print $1}' lineages | sed -n ${SGE_TASK_ID}p)

mkdir raxmlng/${LINEAGE}

module load anaconda3
conda activate raxml-ng

#Convert file for RAxML-NG

raxml-ng --parse \
         --msa alignments/${LINEAGE}/${LINEAGE}_concat.phy \
         --model alignments/${LINEAGE}/${LINEAGE}_partition.txt \
         --prefix raxmlng/${LINEAGE}/${LINEAGE}_concat

#Run ML tree search and bootstrapping for 1000 iterations

raxml-ng --all \
         --msa alignments/${LINEAGE}/${LINEAGE}_concat.phy \
         --model alignments/${LINEAGE}/${LINEAGE}_partition.txt \
         --prefix raxmlng/${LINEAGE}/${LINEAGE}_concat \
         --seed 2 \
         --threads ${NSLOTS} \
         --bs-trees autoMRE{1000}

#Check convergence

raxml-ng --bsconverge \
         --bs-trees raxmlng/${LINEAGE}/${LINEAGE}_concat.raxml.bootstraps \
         --prefix raxmlng/${LINEAGE}/${LINEAGE}_concat_convergence_test \
         --seed 2
