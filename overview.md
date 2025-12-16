# Virome Analysis in Tomato (*Solanum lycopersicum*)using Next Generation Sequencing (NGS) 
*This repository contains scripts and documentation for analyzing virome data in tomato plants using Next Generation Sequencing (NGS) techniques. The analysis focuses on identifying and characterizing viral sequences present in tomato samples.*

### Brief Overview

*This project was completed as a part of my master's degree in Plant Pathology at the University of Agricultural Sciences, Bangalore, India. The main objective was to explore the virome of tomato plants, which are susceptible to various viral infections that can significantly impact crop yield and quality.*

### Introduction

Tomato (*Solanum lycopersicum*) is a widely cultivated vegetable crop that is vulnerable to numerous viral pathogens. Understanding the virome of tomato plants is crucial for developing effective disease management strategies. Conventional diagnostic approaches like ELISA, PCR and Sanger sequencing require prior knowledge of the specific virus being tested, which can hinder the identification of mixed infections and novel viruses. To overcome these limitations, this project utilizes Next Generation Sequencing (NGS) to identify and analyze novel and previously reported viral sequences present in tomato samples collected from southern part of India.

### Data Description

### Methods

Total RNA was extracted from sixty unique tomato leaf samples, collected from different locations using a commercial RNA extraction kit. RNA was pooled from all samples in equimolar concentrations and ribosomal RNA was depleted from the pooled total RNA using the Qiaseq FastSelect Plant rRNA Kit. The RNA-seq library was constructed using the Kapa RNA HyperPlus Kit according to the manufacturer’s protocol. The library was quantified using a Qubit 4.0 Fluorometer with a DNA HS Assay Kit. To determine the insert size of the library, we analyzed it on the TapeStation 4200. Raw data was obtained from sequencing the library using the Illumina NovaSeq 6000 S4 300 cycles kit, generating 2x150 bp paired-end reads.

For this analysis, I will work with paired-end raw sequencing data in FASTQ format obtained from the Illumina High Throughput Sequencing platform. 

Goal: The expected output of this project includes a comprehensive list of identified viral sequences and their taxonomic classification. The results will provide insights into the diversity of viruses present in tomato plants and their potential impact on plant health. Additionally, phylogenetic analysis will help understand the evolutionary relationships among the identified viral sequences.

### Data analysis workflow

**All the files required for the analysis, including raw FASTQ files, scripts, and results, will be found at `/fs/ess/PAS2880/users/sskavya123/FP`**

The analysis workflow consists of the following steps:

1. Determining the quality of raw reads.
2. Trimming of low-quality bases and adapter sequences
3. Mapping reads to the host (tomato) reference genome to filter out host sequences.
4. Extracting unmapped (viral) reads.
5. De novo assembly of filtered viral reads.
6. Identification of viral sequences using BLAST against a viral reference database.
7. Reconstruction of viral genomes.
8. Phylogenetic analysis of identified viral sequences

### Software and Tools Used
- FastQC 

Since the raw data is in FASTQ format, I will use FastQC tool for analysing the quality.
- Trimgalore 

I will use Trimgalore for trimming low-quality bases and adapter sequences from the raw reads. Based on the FastQC reports, I will set appropriate parameters for trimming.
- Star aligner

I will use STAR aligner to map the trimmed reads to the tomato reference genome to filter out host sequences.
- SPAdes assembler 

Since I am aiming to identify novel viral sequences, I will use SPAdes assembler for de novo assembly of the filtered viral reads.
- BLAST 

I will use BLAST to compare the assembled sequences against the viral reference database.
- Clustal Omega

Since the output of *de novo* assembly is a set of contigs, I will use Clustal Omega for multiple sequence alignment of the identified viral sequences.
- Bioedit 

I will use Bioedit for manual editing and reconstruction of viral genomes from the aligned sequences.
- iqtree 

I will use iqtree for phylogenetic analysis of the identified viral sequences.

**This is a work in progress, and I am still learning to use some of these tools effectively. I am also exploring alternative tools that may better suit certain steps of the analysis.**

*Singularity containers will be created for each software tool to ensure reproducibility and ease of deployment across different computing environments.*

### Structure of the analysis

I will use a protocol markdown file to document the entire analysis. It will include consequtive steps, commands used, parameters set, and explanations for each step. 

### Directory structure

I will create a data directory to store the raw FASTQ files. I will write scripts for each step of the analysis and save them in a scripts directory. A dedicated results directory will be created to store the output files generated at each step of the analysis.

*The data directory* -  `/fs/ess/PAS2880/users/sskavya123/FP/data` will be a read -only directory to ensure the integrity of the raw data. The scripts and results directories will have appropriate permissions to allow for modifications and updates as needed.

*The scripts directory*-  `/fs/ess/PAS2880/users/sskavya123/FP/scripts` will contain individual scripts - `fastqc.sh, trimgalore.sh, star_index.sh, star_align.sh, spades.sh, blast.sh, clustalo.sh and iqtree.sh`. I will run the scripts as batch jobs and the relevant sbatch command will be in the protocol markdown file. Each script will take input files and output files as arguments to ensure flexibility and reusability. The scripts will also contain a provision for saving log files to track the progress and any errors encountered during execution.

*The results directory*-  `/fs/ess/PAS2880/users/sskavya123/FP/results` will be organized in subdirectories within the results directory, corresponding to each step of the analysis. For example, there will be separate subdirectories for FastQC reports, trimmed reads, assembled contigs, BLAST results, and phylogenetic trees.

### Aspects of project still uncertain about:

We have learnt FastQC, Trimgalore, STAR aligner in detail. However, I am still learning to use SPAdes assembler, BLAST, Clustal Omega and iqtree. I am concerned about leveraging the output of BLAST to Clustal Omega for multiple sequence alignment since I might have to extract specific sequences based on best hits from BLAST results. I am also concerned about ambiguous bases in the reconstructed viral genomes and how to handle them effectively.

I am hoping to finish the project. This will give me a strong foundation in NGS data analysis and virome research, which will be valuable for my future career in genomics.

### Future Directions

I am planning to create publication quality figures and tables summarizing the findings of the virome analysis using R. I am also considering exploring additional analyses, such as functional annotation of viral genes and variant analysis of viral sequences to identify potential mutations associated with virulence or host adaptation.

### Reasons for choosing this project

I chose this data set because I generated it myself during master's degree and I still could not publish it because I was not confident enough with the data analysis part. I hope to publish this work after completing this project.