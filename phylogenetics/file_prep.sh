cat GenePull/${LINEAGE}/*${MARKER}.fa alignments/${LINEAGE}/${MARKER}.fasta > tmp && mv tmp alignments/${LINEAGE}/${MARKER}.fasta

sed 's/>[^ ]* />/' GAPDH.fasta | sed 's/glyce.*//' | sed 's/GAPDH.*//' > tmp && mv tmp GAPDH.fasta
sed 's/>[^ ]* />/' HIS3.fasta | sed 's/hist.*//' | sed 's/HIS.*//' > tmp && mv tmp HIS3.fasta
sed 's/>[^ ]* />/' TUB2.fasta | sed 's/beta.*//' > tmp && mv tmp TUB2.fasta

sed 's/ /_/g' GAPDH.fasta | sed 's/_$//' | sed 's/:/_/' | sed 's/(//' | sed 's/)//' > tmp && mv tmp GAPDH.fasta
