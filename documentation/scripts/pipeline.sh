#!/bin/bash

# Check if the correct number of arguments is provided
if [ $# -ne 2 ]; then
  echo "Usage: $0 <folder_path> <program_folder>"
  exit 1
fi

# Assign the input arguments to variables
folder_path=$1
program_folder=$2
ref_genome="${program_folder}/reference/ref-genome/GRCh38.p13.genome.fa"

# Check if the folder exists
if [ ! -d "$folder_path" ]; then
  echo "Folder does not exist: $folder_path"
  exit 1
fi

# Check if the program folder exists
if [ ! -f "$program_folder" ]; then
  echo "Program folder does not exist: $program_folder"
  echo "Download the program_folder from https://github.com/rpianezza/Human-Genes"
  exit 1
fi

# Create a "mapped" subfolder if it doesn't exist
mapped_folder="$folder_path/mapped"
if [ ! -d "$mapped_folder" ]; then
  mkdir "$mapped_folder"
fi

# Create a "output" subfolder if it doesn't exist
output_folder="$folder_path/output"
if [ ! -d "$output_folder" ]; then
  mkdir "$output_folder"
fi

# Loop through the fastq files in the folder
for fastqgz_file in "$folder_path"/*.fastq.gz; do
  # Check if any fastq files exist
  if [ ! -e "$fastqgz_file" ]; then
    echo "No fastq.gz files found in $folder_path"
    exit 1
  fi

  # Extract the file name without the extension
  filename=$(basename "$fastqgz_file" .fastq.gz)

  # Run BWA-MEM to map fastq to fasta
  gunzip -k "$fastqgz_file"
  bwa mem "$ref_genome" "${folder_path}/${filename}.fastq" > "${mapped_folder}/${filename}.sam"
  samtools view -bS -o "${mapped_folder}/${filename}.bam" "${mapped_folder}/${filename}.sam" > "/dev/null"
  samtools sort "${mapped_folder}/${filename}.bam" -o "${mapped_folder}/${filename}.sorted.bam"
  samtools index "${mapped_folder}/${filename}.sorted.bam"
  echo "Fastq file mapped"

  # Call script to calculate SCGs coverage
  bash "${program_folder}/scripts/scg-coverage.sh "${mapped_folder}/${filename}.sorted.bam" "${program_folder}/annotation/genes/scg.bed" "${output_folder}/${filename}.scg.tsv"
  

  rm "${folder_path}/${filename}.fastq"
  rm "${mapped_folder}/${filename}.sam"
  rm "${mapped_folder}/${filename}.bam"
done


