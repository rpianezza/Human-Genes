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
if [ ! -d "$program_folder" ]; then
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

  # Run BWA-MEM to map fastq to bam
  echo "Extracting ${filename}.fastq.gz..."
  gunzip -k "$fastqgz_file"
  echo "Mapping ${filename}.fastq... (this will take a while...)"
  bwa mem "$ref_genome" "${folder_path}/${filename}.fastq" > "${mapped_folder}/${filename}.sam" 2> "/dev/null"
  samtools view -q 10 -bS -o "${mapped_folder}/${filename}.bam" "${mapped_folder}/${filename}.sam" > "/dev/null"
  samtools sort "${mapped_folder}/${filename}.bam" -o "${mapped_folder}/${filename}.sorted.bam"
  samtools index "${mapped_folder}/${filename}.sorted.bam"

  # Call scripts to calculate SCGs coverage and PCGs exons coverage
  echo "Calculating SCGs mean coverage..."
  bash "${program_folder}/scripts/scg-coverage.sh" "${mapped_folder}/${filename}.sorted.bam" "${program_folder}/annotation/genes/scg.bed" "${output_folder}/${filename}.scg.tsv"
  echo "Calculating PCGs exons coverage... (this will take a while...)"
  bash "${program_folder}/scripts/pcg-coverage.sh" "${mapped_folder}/${filename}.sorted.bam" "${program_folder}/annotation/exons/merged-exons-pcg.bed" "${output_folder}/${filename}.pcg-exons.tsv"
  
  # Merge PCG exons coverage
  echo "Merging PCGs exons..."
  python "${program_folder}/scripts/pcg-merged-coverage.py" "${output_folder}/${filename}.pcg-exons.tsv" "${output_folder}/${filename}.pcg.tsv"
  
  # Normalize PCG coverage using mean SCG coverage
  echo "Normalizing PCGs copynumber..."
  python "${program_folder}/scripts/normalize.py" "${output_folder}/${filename}.scg.tsv" "${output_folder}/${filename}.pcg.tsv" "${output_folder}/${filename}.pcg-normalized.tsv"
  
  # Normalize GC content
  echo "Correcting for GC-bias..."
  mkdir "${output_folder}/GC-plots"
  Rscript "${program_folder}/scripts/gc-normalization.R" --gc_scg "${program_folder}/gc-content/scg-gc.tsv" --cov_scg "${output_folder}/${filename}.scg.tsv" --pcg_cn "${output_folder}/${filename}.pcg-normalized.tsv"  --gc_pcg "${program_folder}/gc-content/pcg-gc.tsv" --output "${output_folder}/"

  rm "${folder_path}/${filename}.fastq"
  rm "${mapped_folder}/${filename}.sam"
  rm "${mapped_folder}/${filename}.bam"
  rm "${output_folder}/${filename}.scg.tsv"
  rm "${output_folder}/${filename}.pcg.tsv"
  rm "${output_folder}/${filename}.pcg-normalized.tsv"
done