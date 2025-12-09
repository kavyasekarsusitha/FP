#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-blast-%j.out
set -euo pipefail

#constants
BLAST=oras://community.wave.seqera.io/library/blast:2.17.0--3d1eb1104ccfd59c

# Copy the placeholder variables
QUERY="$1"
REF_DB="$2"
OUTDIR="$3"

# Initial logging
echo "# Starting script blast.sh"
date
echo "# Input FASTQ file:   $QUERY"
echo "# Reference DB:      $REF_DB"
echo "# Output dir:         $OUTDIR"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$OUTDIR"/logs

# Run BLAST
echo "Running BLAST..."
apptainer exec "$BLAST" blastn -query "$QUERY" \
       -db "$REF_DB" \
       -max_target_seqs 1 \
       -outfmt "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send slen stitle" \
       -evalue 1e-5 \
       -num_threads 8 \
       > "$OUTDIR/blastn.outfmt6"

# Final logging
echo "# Finished script blast.sh"
date