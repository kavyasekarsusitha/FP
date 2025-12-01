#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-star_align-%j.out
set -euo pipefail

# Constants
STAR_CONTAINER=oras://community.wave.seqera.io/library/star:2.7.11b--84fcc19fdfab53a4

# Copy the placeholder variables
R1="$1"
R2="$2"
index_dir="$3"
gtf="$4"
outdir="$5"

# Infer the "sample ID" - we need this for the output file specification
sample_id=$(basename "$R1" _R1.fastq.gz)

# Report start of script and variables
echo "# Starting script star_align.sh"
date
echo "# Input R1 FASTQ file:            $R1"
echo "# Input R2 FASTQ file:            $R2"
echo "# Input genome index:             $index_dir"
echo "# Input annotation GTF:           $gtf"
echo "# Output dir:                     $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run STAR
apptainer exec "$STAR_CONTAINER" STAR \
    --readFilesIn "$R1" "$R2" \
    --genomeDir "$index_dir" \
    --sjdbGTFfile "$gtf" \
    --alignIntronMin 40 \
    --alignIntronMax 23000 \
    --outFilterMultimapNmax 10 \
    --runThreadN 8 \
    --readFilesCommand zcat \
    --outFileNamePrefix "$outdir"/"$sample_id"_ \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx

# --alignIntronMin 40 option is used to set the minimum intron length to 40 bases.
# --alignIntronMax was set to 23000 to accommodate larger introns.
# --outFilterMultimapNmax was set to default of 10 to allow reads to map to up to 10 locations.
# --readFilesCommand zcat                => Tell STAR that the FASTQ files are gzipped
# --outSAMtype BAM SortedByCoordinate    => Request a sorted BAM file as output (instead of unsorted SAM)
# --outFileNamePrefix "$outdir"/"$sample_id"_ => Specify not just an outdir but a "sample ID" prefix
# --outReadsUnmapped Fastx               => Output unmapped reads in FASTQ format for downstream viral analysis

# Final logging
echo
echo "# Used STAR version:"
apptainer exec "$STAR_CONTAINER" STAR --version
echo "# Successfully finished script star_align.sh"
date
