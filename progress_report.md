# PROGRESS REPORT
### Virome Analysis in Tomato (*Solanum lycopersicum*)using Next Generation Sequencing (NGS) 
*This report summarizes the progress made in the virome analysis project focusing on tomato plants using NGS techniques. The analysis aims to identify and characterize viral sequences present in tomato samples.*

The protocol_README.md file contains the detailed steps followed for the analysis. README file submitted as proposal last week contains the rationale behind using specific tools and methods for the analysis.

### Work Completed So Far
1. I have completed running fastqc on the raw FASTQ files to assess the quality of the sequencing data. I have uploaded the per base sequence quality plots in the github repository. I am concerned about the quality drop towards the end of the reads, which might affect downstream analyses.

The script fastqc.sh can be fpound in the scripts/ directory. I have run the script and I am confident that it has executed successfully without errors.

2. Up next, I completed trimming the adapters and low-quality bases using Trimgalore. The trimmed reads have been saved in the results/trimgalore directory. I have checked the quality of the trimmed reads using FastQC again, and the quality has improved significantly.

Since I have one pair of FASTQ files, I did not run MultiQC for aggregating the FastQC reports. I am confident that the trimming step was successful.

3. Then I planned to filter the virus associated reads before proceeding to *de novo* assembly. 

I downloaded the .fna and .gtf files for the tomato reference genome from NCBI and moved them to the data/ref directory in my working directory. Then I used star for creating the genome index - correspoding script star_index.sh can be found in the scripts/ directory. I ran the script and it executed successfully.

Some options I used for STAR indexing which are different from default values:

since the sequencing read length is 150, sjdbOverhang is set to 149
--genomeSAindexNbases 13 is set considering the genome size for fine-grained index

The script was executed successfully without errors.

4. Next, I aligned the trimmed reads to the tomato reference genome using STAR aligner to filter out host reads. The script star_align.sh can be found in the scripts/ directory.

Some options I used for STAR alignment which are different from default values:

--alignIntronMin 40 : I found a research article  `doi: 10.1186/s12864-021-08212-x` which stated that the mean intron length was around 632 bp for tomato genome. But, I could not find the exact range. I later found a conversation thread about STAR usage `https://groups.google.com/g/rna-star/c/R2A352vfhlI` and found a suggestion to set the minimum intron length to 40 for tomato genome. Hence, I set --alignIntronMin to 40.

--alignIntronMax 23000 : Similarly for maximum intron length, I set 23000 as suggested in the same thread.

--outFilterMultimapNmax 10 : To allow for reads mapping to multiple locations up to 10, seemed reasonable considering the complexity of the tomato genome.

--outReadsUnmapped Fastx : To save the unmapped reads in FASTQ format for downstream virome analysis.

The alignment step was completed successfully without errors. The unmapped reads have been saved in the data/unmapped directory.

5. I used SPAdes for *de novo* assembly of the unmapped reads to reconstruct potential viral genomes. The script spades_assembly.sh can be found in the scripts/ directory. 

Some options I used for SPAdes assembly which are different from default values:

--only-assembler tells SPAdes to skip the read error correction step and proceed directly to assembly. Since the reads have already been quality trimmed using TrimGalore, I opted to skip error correction to save computational time.

The script was executed successfully without errors. The assembled contigs have been saved in the results/spades directory.

To-do List:

1. Next, I plan to annotate the assembled contigs using BLASTn against the NCBI viral database to identify potential viral sequences. I will create a BLAST database from the viral sequences and run BLASTn to find matches for the assembled contigs. 

2. An alternative approach is to consider using DIAMOND for protein-level annotation, which might help identify more divergent viral sequences that may not have close nucleotide matches in the database. But I will first try BLASTn and see how many viral contigs I can identify.

3. Following annotation, I will reconstruct viral genomes from the identified viral contigs using Clustal Omega and BioEdit for multiple sequence alignment and manual editing, respectively.

4. Finally, I will compile the results, including identified viral sequences, their annotations, and any relevant statistics, into a comprehensive report.

Feedback Request:

I would appreciate feedback on the following aspects of my analysis so far:

1. I have changed some default parameters in STAR alignment based on literature and forum suggestions. Are these changes appropriate for my analysis, or should I consider alternative settings?

2. I skimmed the spades contigs.fasta file and found some host sequences. Should I consider additional filtering steps before assembly to further reduce host contamination? Please refer to T1_Log.final.out file in results/star directory for alignment statistics. 
