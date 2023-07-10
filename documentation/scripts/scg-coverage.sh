#!/bin/bash

#This script takes as input a BED file of exons with these columns: "chr", "start", "end", "gene".
#It calculates the mean coverage of each exons and output a tsv file ("exon", "coverage", "gene").

# Check if the required tools are installed
command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools is required but not installed. Aborting."; exit 1; }

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_file> <bed_file> <output_file>"
    exit 1
fi

bam_file=$1
bed_file=$2
output_file=$3

# Iterate over each entry in the BED file
while read -r bed_entry; do
    # Extract chromosome, start, end positions, and gene from the BED entry
    chromosome=$(echo "$bed_entry" | awk '{print $1}')
    start=$(echo "$bed_entry" | awk '{print $2}')
    end=$(echo "$bed_entry" | awk '{print $3}')
    gene=$(echo "$bed_entry" | awk '{print $4}')

    # Calculate mean coverage for gene using samtools depth
    mean_coverage=$(samtools depth -r "$chromosome:$start-$end" "$bam_file" | awk '{sum += $3} END {if (NR > 0) print sum / NR; else print 0}')

    # Print the result for scg
    printf "%s\t%s\t%s\t%s\n" "$chromosome:$start-$end" "$mean_coverage" >> "$output_file"

done < "$bed_file"