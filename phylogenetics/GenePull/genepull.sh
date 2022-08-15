#!/bin/sh

module bedtools
module blast+

 ~/Scripts/GenePull -g HIS3_example.fasta -a ../../../assessment/366226/blobtools/366226_spades_polished_filtered_nocontam.fa -o 366226_HIS3

sed -i 's/>.*$/>355084/' 355084_GAPDH.fa
