#!/bin/bash


TRANSCRIPTOME="gencode.v47.transcripts.fa.gz"    # Path to transcriptome FASTA file
INDEX_DIR="salmon_index"                         # Directory where the index will be created
DECOYS_FILE="decoys.txt"                         # Path to decoys.txt file
KMER_SIZE=31                                     # K-mer size


# build the Salmon index without decoy
salmon index -t $TRANSCRIPTOME -i $INDEX_DIR --gencode -p 15 -k $KMER_SIZE

echo "Indexing completed."
