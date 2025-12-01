#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-spades_assembly-%j.out
set -euo pipefail

# Constants
SPADES_CONTAINER=oras://community.wave.seqera.io/library/spades:4.2.0--3313822b80929818

# Copy the placeholder variables
R1="$1"
R2="$2"
outdir="$3"


# Report start of script and variables
echo "# Starting script spades.sh"
date
echo "# Input R1 FASTQ file:            $R1"
echo "# Input R2 FASTQ file:            $R2"
echo "# Output dir:                     $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run SPAdes
apptainer exec "$SPADES_CONTAINER" spades.py \
    --pe1-1 "$R1" \
    --pe1-2 "$R2" \
    -o "$outdir"/spades \
    -t 8 \
    --only-assembler
    
# --pe1-1 and --pe1-2 specify the paired-end reads.
# -o specifies the output directory for the assembly results.
# -t 8 specifies the number of threads to use for the assembly.
# --only-assembler tells SPAdes to skip the read error correction step and proceed directly to assembly.
# Final logging
echo
echo "# Used SPAdes version:"
apptainer exec "$SPADES_CONTAINER" spades.py --version
echo "# Successfully finished script spades.sh"     
