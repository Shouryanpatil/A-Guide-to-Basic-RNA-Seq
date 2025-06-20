#!/bin/bash

# Create output directories if they don't exist
mkdir -p aligned_bam/raw_bam
mkdir -p aligned_bam/sorted_bam

# Define sample names
samples=(
  Sample_1_1
  Sample_1_2
  Sample_1_3
  Sample_1_4
  Sample_2_1
  Sample_2_2
  Sample_2_3
  Sample_2_4
)

# Loop through each sample and convert, sort, index
for sample in "${samples[@]}"; do
  echo "Processing $sample..."

  # Convert SAM to BAM
  samtools view -S -b aligned_bam/aligned_sam/${sample}.sam > aligned_bam/raw_bam/${sample}.bam

  # Sort BAM
  samtools sort aligned_bam/raw_bam/${sample}.bam -o aligned_bam/sorted_bam/${sample}_sorted.bam

  # Index BAM
  samtools index aligned_bam/sorted_bam/${sample}_sorted.bam

  echo "$sample done."
  echo "-----------------------------"
done

echo "All samples processed."
