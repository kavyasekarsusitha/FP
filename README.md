# Virome Analysis in Tomato (*Solanum lycopersicum*)using Next Generation Sequencing (NGS) 
### Brief Overview

*This project was completed as a part of my master's degree in Plant Pathology at the University of Agricultural Sciences, Bangalore, India. The main objective was to explore the virome of tomato plants, which are susceptible to various viral infections that can significantly impact crop yield and quality.*

### Introduction

Tomato (*Solanum lycopersicum*) is a widely cultivated vegetable crop that is vulnerable to numerous viral pathogens. Understanding the virome of tomato plants is crucial for developing effective disease management strategies. Conventional diagnostic approaches like ELISA, PCR and Sanger sequencing require prior knowledge of the specific virus being tested, which can hinder the identification of mixed infections and novel viruses. To overcome these limitations, this project utilizes Next Generation Sequencing (NGS) to identify and analyze novel and previously reported viral sequences present in tomato samples collected from southern part of India.

### Methods

#### Data analysis workflow
**All the files required for the analysis, including raw FASTQ files, scripts, and results, will be found at `/fs/ess/PAS2880/users/sskavya123/FP`**

**The protocol_README.md file contains detailed instructions on how to run the analysis pipeline. The scripts used in the analysis are available in the `scripts/` directory.**


The analysis workflow consists of the following steps:
1. Determining the quality of raw reads - FastQC
2. Trimming of low-quality bases and adapter sequences - Trimgalore
3. Mapping reads to the host (tomato) reference genome to filter out host sequences - Star aligner
4. Extracting unmapped (viral) reads 
5. De novo assembly of filtered viral reads - SPAdes assembler
6. clustering similar contigs to reduce redundancy - CD-HIT
7. Creating a local BLAST database of viral genomes - makeblastdb
8. Identification of viral sequences using BLAST against a viral reference database - BLASTn
9. Downstream analysis - R

### Results

The results of the analysis include a comprehensive list of identified viral sequences along with their taxonomic classification. The findings provide insights into the diversity of viruses present in tomato plants and their potential impact on plant health. 

Virome analysis identified 91 unique viral sequences, containing 76 RNA viruses, 12 DNA viruses, and 3 satellite viruses. Notably, several economically significant viruses were detected, including Tomato yellow leaf curl virus (TYLCV), Chilli vienal mottle virus, and Cucumber mosaic virus (CMV).

