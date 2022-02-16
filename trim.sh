#!/bin/bash
mkdir qc_reads
SAMPLES=(
  "ACBarrie"
  "Alsen"
  "Baxter"
  "Chara"
  "Drysdale"
)

for SAMPLE in "${SAMPLES[@]}"; do

  fastp \
    -i "$PWD/${SAMPLE}_R1.fastq.gz" \
    -I "$PWD/${SAMPLE}_R2.fastq.gz" \
    -o "qc_reads/${SAMPLE}_R1.fastq.gz" \
    -O "qc_reads/${SAMPLE}_R2.fastq.gz" \
    --html "qc_reads/${SAMPLE}_fastp.html" 
done
