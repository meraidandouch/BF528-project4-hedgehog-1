#!/bin/bash
GENOME_FA="GRCh38.primary_assembly.genome.fa.gz" # path to genome FASTA file 

echo "Using genome file: $GENOME_FA"
# decoys.txt from the genome FASTA to obtain the cleaned chromosome names in the file  
zgrep "^>" $GENOME_FA | sed 's/>//;s/ .*//' > decoys.txt
echo "Decoy extraction completed. Output saved to decoys.txt"
