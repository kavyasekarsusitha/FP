#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-cd_hit-%j.out
set -euo pipefail

# Constants
CDHIT_CONTAINER=oras://community.wave.seqera.io/library/cd-hit:4.8.1--b8ee5a3722ee5d85

# Copy the placeholder variables

input="$1"
outdir="$2"


# Report start of script and variables
echo "# Starting script cd_hit.sh"
date
echo "# Input FASTA file:               $input"
echo "# Output dir:                     $outdir"


# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run CD-HIT
apptainer exec "$CDHIT_CONTAINER" cd-hit \
    -i "$input" \
    -o "$outdir"/cdhit_output.fasta \
    -c 0.90 \
    -n 5 \
    -T 8 \
    -M 0  \
    -d 0     
