#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p aligned_bam/aligned_sam

# Define genome index path
GENOME_INDEX="genome_index/grch38/genome"

# Sample name array
samples=(
  "Sample_1_1"
  "Sample_1_2"
  "Sample_1_3"
  "Sample_1_4"
  "Sample_2_1"
  "Sample_2_2"
  "Sample_2_3"
  "Sample_2_4"
)

# Loop through each sample and run HISAT2
for sample in "${samples[@]}"
do
  echo "Aligning $sample..."
  hisat2 -p 4 -x $GENOME_INDEX \
    -U trimmed_fastq/${sample}_trimmed.fastq \
    -S aligned_bam/aligned_sam/${sample}.sam
  echo "$sample alignment complete."
done

echo "All alignments completed successfully."
