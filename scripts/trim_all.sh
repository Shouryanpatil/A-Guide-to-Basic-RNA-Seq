#!/bin/bash

# Make sure output directory exists
mkdir -p trimmed_fastq

# Array of input files and corresponding sample names
declare -a samples=(
  "SRR33341769 Sample_1_1"
  "SRR33341768 Sample_1_2"
  "SRR33341767 Sample_1_3"
  "SRR33341766 Sample_1_4"
  "SRR33341765 Sample_2_1"
  "SRR33341764 Sample_2_2"
  "SRR33341763 Sample_2_3"
  "SRR33341762 Sample_2_4"
)

# Path to adapter file (adjust if needed)
ADAPTER="TruSeq3-SE.fa"

# Loop through each sample and run Trimmomatic
for entry in "${samples[@]}"; do
  read SRR_ID SAMPLE_NAME <<< "$entry"
  
  echo "Trimming $SAMPLE_NAME ($SRR_ID)..."
  
  trimmomatic SE -threads 4 -phred33 \
    "${SRR_ID}.fastq" \
    "trimmed_fastq/${SAMPLE_NAME}_trimmed.fastq" \
    ILLUMINACLIP:$ADAPTER:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
  
  echo "Finished $SAMPLE_NAME"
done

echo "All samples trimmed."
