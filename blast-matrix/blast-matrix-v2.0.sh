#!/bin/bash

# Set the paths and filenames
BLAST_PATH="/usr/local/Caskroom/miniconda/base/bin/blastn"
QUERY_FILE="$1"
SUBJECT_FILE="$2"
OUTPUT_FILE="$3"

# Read the query sequence from the input FASTA file
query_id=$(grep -oE "^>.+" "$QUERY_FILE" | tr -d ">")
query_sequence=$(grep -v ">" "$QUERY_FILE")

# Perform the BLAST search and save the output
echo -e ">$query_id\n$query_sequence" | $BLAST_PATH -query - -subject <(cat "$SUBJECT_FILE") -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" >> "$OUTPUT_FILE"

echo "BLAST of $query_id completed"