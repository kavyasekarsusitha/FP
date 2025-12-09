# Virome Analysis in Tomato (*Solanum lycopersicum*)using Next Generation Sequencing (NGS) 

- Author: Kavya Sekar Susitha
- Date: 26 November, 2025
- GitHub repository: 
- This was run at the pitzer cluster of the Ohio Supercomputer Center (<www.osc.edu>)
- OSC working dir: /fs/ess/PAS2880/users/sskavya123/FP

## Prerequisites

- **Input files**:
  - The input files are:
    - Paired-end FASTQ files T1_R1.fastq.gz and T1_R2.fastq.gz
    - Tomato reference genomes files (assembly and annotation)
- **Output Directories**:
  - The scripts will create output directories as needed 

- **Software**:
  - The scripts will run programs using Singularity/Apptainer containers
    available online, with links contained inside the scripts.

## Set up: Get the input files

1. I have downloaded the tomato reference genome assembly and annotation files from NCBI and uploaded it in my directory `/fs/ess/PAS0471/kavya/GCF_036512215.1`. I will move the files to my working directory within data/ref subdirectory.

   ```bash
   mkdir -p data/ref
   mv /fs/ess/PAS0471/kavya/GCF_036512215.1 data/ref/
   ```

## Set up: Set file paths

```bash
# Inputs:
fastq_dir=data/fastq
ref_assembly=data/ref/GCF_036512215.1/GCF_036512215.1_SLM_r2.1_genomic.fna
ref_annotation=data/ref/GCF_036512215.1/genomic.gtf

# Base output dir:
outdir=results
```

## Step 1: Read QC with FastQC

The `fastqc.sh` script runs on 1 FASTQ file at a time and takes 2 arguments -
the input FASTQ file and output dir:

```bash
for fastq in "$fastq_dir"/*fastq.gz; do
    sbatch scripts/fastqc.sh "$fastq" "$outdir"/fastqc
done
```
After checking the log files and confirming successful completion, I will move the slurm log files to the logs directory within results/fastqc.

```bash
mv -v slurm-fastqc* results/fastqc/logs/
```

Since I have only two fastq files, I will not use multiqc for aggregating the fastqc reports.

## Step 2: Read trimming with TrimGalore

The script `trimgalore.sh` takes 3 arguments - input R1 FASTQ file, input R2 FASTQ file, and output dir. 

I will set the default quality threshold to 20 using `--quality 20` parameter and adapter trimming will be done automatically by TrimGalore.

```bash
for R1 in "$fastq_dir"/*_R1.fastq.gz; do
    R2="${R1/_R1/_R2}"
    sbatch scripts/trimgalore.sh "$R1" "$R2" "$outdir"/trimgalore
done
```
moving the slurm log files to the logs directory within results/trimgalore.

```bash
mv slurm-trimgalore* results/trimgalore/logs
```
## Step 3: Creating STAR genome index
The `star_index.sh` script creates the STAR genome index and takes 3 arguments - reference genome assembly file, reference annotation file, and output dir for the index.

```bash
sbatch scripts/star_index.sh "$ref_assembly" "$ref_annotation" results/star/index
```
After confirming successful completion, I will move the slurm log files to the logs  directory within results/star/index.

```bash
mv slurm-star_index* results/star/index/logs
``` 
## Step 4: Aligning reads to the host genome
The `star_align.sh` script aligns the trimmed reads to the host genome and takes 4 arguments - input R1 FASTQ file, input R2 FASTQ file, STAR index dir, and output dir.

```bash
for R1 in results/trimgalore/*R1.fastq.gz; do
    R2="${R1/_R1/_R2}"
    sbatch scripts/star_align.sh "$R1" "$R2" results/star/index "$ref_annotation" results/star
done
```
After confirming successful completion, I will move the slurm log files to the logs directory within results/star.

```bash
mv slurm-star_align* results/star/logs
```
## Step 5: *de novo* assembly of unmapped reads

The `spades.sh` script assembles the unmapped reads using SPAdes and takes 3 arguments - input R1 FASTQ file, input R2 FASTQ file, and output dir.

```bash
sbatch scripts/spades.sh data/unmapped/T1_R1_unmapped.fastq data/unmapped/T1_R2_unmapped.fastq results/spades
```
After confirming successful completion, I will move the slurm log files to the logs directory within results/spades.

```bash
mv slurm-spades* results/spades/logs
```
## Step 6: cd-hit clustering of contigs
The `cdhit.sh` script clusters the contigs using cd-hit and takes 2 arguments - input contigs file and output dir.

```bash
sbatch scripts/cd_hit.sh results/spades/spades/contigs.fasta results/cdhit
```
After confirming successful completion, I will move the slurm log files to the logs directory within results

```bash
mv slurm-cd_hit* results/cdhit/logs
```

## Step 7: Create BLAST database 
The `blast_db.sh` script creates a BLAST database from viral reference genomes.

As the command did not have provisions for mentioning output directory, I have modified the script to create and move to the desired output directory before downloading the viral reference genomes.
```bash
sbatch scripts/blast_db.sh
```
```bash
mv slurm-blast_db* results/db/viral_refseq/logs
``` 

## step 8: Run BLASTn of contigs against viral reference database

# The `blastn.sh` script runs BLASTn of the clustered contigs against the viral reference database and takes 3 arguments - input contigs file, reference database, and output dir.

```bash
sbatch scripts/blast.sh results/cdhit/cdhit_output.fasta results/db/viral_refseq/viral_refseq_db results/blastn
``` 
After confirming successful completion, I will move the slurm log files to the logs directory within results/blastn.

```bash
mv slurm-blast* results/blastn/logs
```