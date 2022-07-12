# Endophyte Genomes
 
![Pipeline workflow](pipeline.png)

The pipeline was written for and run on Queen Mary University of London's [Apocrita HPC facility](http://doi.org/10.5281/zenodo.438045) which uses the Univa Grid Engine batch-queue system. This means that many of the bash scripts (`.sh` file endings) specify core allocation, run times and memory usage allocation that may need to be adapted for different platforms.

---

## 1 Read trimming/basecalling

`cd reads`

### Short-reads

1. `qsub trimmomatic.sh` trims raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic); requires `TruSeq3-PE.fa` file with adapter sequences downloaded from [here](https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa) (for Illumina NovaSeq 6000 151bp paired-end reads).
2. `qsub fastqc.sh` checks trimmed read quality with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### Long-reads

`qsub guppy.sh` performs fast basecalling of raw MinION read data.

## 2 *De novo* genome assembly

`cd denovo_assembly`

1. `./submit_assembly.sh` makes new directory and submits job scripts for each assembly tool - short-read tools `abyss.sh` ([ABySS](https://github.com/bcgsc/abyss)), `megahit.sh` ([MEGAHIT](https://github.com/voutcn/megahit)) and `spades.sh` ([SPAdes](https://github.com/ablab/spades)), and long-read tools `flye.sh` ([Flye]https://github.com/fenderglass/Flye), `raven.sh` ([Raven]https://github.com/lbcb-sci/raven) and `spades_hybrid` (hybridSPAdes). For all tools except for ABySS, these scripts also include short-read mapping with [BWA-MEM](https://github.com/lh3/bwa) for polishing with [Pilon](https://github.com/broadinstitute/pilon) and remove sequences <200bp using [Seqtk](https://github.com/lh3/seqtk) for NCBI compliance.
2. `qsub abyss_comp.sh` compares the assembly stats to choose 'best' kmer size for ABySS (must be done after `abyss.sh` has finished for all kmer sizes and strains), followed by short-read polishing with Pilon.


## 3 Assessment

`cd assessment`

### Assembly tool comparison

1. `./assessment_submit` submits scripts for assembly quality statistics - `quast.sh` ([QUAST](https://github.com/ablab/quast)) and `busco.sh` ([BUSCO](https://busco.ezlab.org/)), which requires the ascomycota_odb10.2020-09-10 BUSCO dataset downloaded from [here](https://busco-data.ezlab.org/v4/data/lineages/)) - and `blast.sh` ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) and `read_mapping.sh`, which maps reads with BWA-MEM, to produce input for BlobTools.
2. `qsub blobtools.sh` submits `blobtools.sh` to run [BlobTools](https://github.com/DRL/blobtools) (must be done after `blast.sh` and `read_mapping.sh` have finished for the strain(s) in question).

### Contamination filtering

1. `qsub remove_contam.sh` removes contigs which BlobTools flagged as belonging to the wrong taxonomic class.
2. `qsub read_mapping_final.sh` - performs a final round of read mapping and produces mapping statisics with [SAMtools](http://www.htslib.org/) to calculate both short- and long-read coverage.

### Final quality statistics

1. `qsub quast_final.sh` reruns QUAST on the contaminant-filtered assemblies.
2. `qsub busco_final.sh` reruns BUSCO on the contaminant-filtered assemblies.
