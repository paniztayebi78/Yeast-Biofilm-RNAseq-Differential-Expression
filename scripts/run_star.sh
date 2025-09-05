≈#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

module load star

FASTQ_DIR="/scratch/ptayebi/BINF6110_Project2/raw_data"
GENOME_DIR="/scratch/lukens/Assignment_2_Genome"
OUTPUT_DIR="/scratch/ptayebi/BINF6110_Project2/aligned_data"

mkdir -p "$OUTPUT_DIR"

for FILE in ${FASTQ_DIR}/*.fastq*; do
    # Get base name (e.g., "SRR10551657")
    base=$(basename "$FILE" | sed 's/.fastq.*//')

    # Set decompression command
    if [[ $FILE == *.gz ]]; then
        CMD="zcat"
    else
        CMD="cat"
    fi

    # Run STAR (SINGLE-END mode)
    STAR --runThreadN 8 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$FILE" \
        --readFilesCommand "$CMD" \
        --outFileNamePrefix "${OUTPUT_DIR}/${base}_" \
        --outSAMtype BAM SortedByCoordinate
done
