#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-blast_db-%j.out
set -euo pipefail

# -------------------------
# Simple viral reference download + BLAST pipeline
# -------------------------

# constants
BLAST=oras://community.wave.seqera.io/library/blast:2.17.0--3d1eb1104ccfd59c

mkdir -p results/db/viral_refseq
cd results/db/viral_refseq

# Download RefSeq viral genomes
echo "Downloading viral reference genomes..."
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.*.genomic.fna.gz

# Unzip
echo "Unzipping..."
gunzip -f viral.*.genomic.fna.gz

# Combine sequences
echo "Combining reference FASTAs..."
cat viral.*.genomic.fna > viral_refseq.fna

# Make BLAST database
echo "Creating BLAST database..."
apptainer exec $BLAST makeblastdb -in viral_refseq.fna -dbtype nucl -out viral_refseq_db

echo "Done!"

