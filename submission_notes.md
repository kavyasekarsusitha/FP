## Instructions to course teachers

1. This file contains notes regarding the submission of the final project. It is intended for course teachers Dr Jelmer and Dr Menuka to understand any specific details about the submission.

2. The protocol_README.md file contains the detailed steps followed in the analysis, along with the relevant code snippets. Each step includes the commands used to execute the analysis, as well as instructions for organizing the output files and logs.

3. The raw data files including the annotations were already downloaded manually and uploaded in the OSC student account prior to starting the course, so that step is not included in the protocol_README.md file. I simply moved the files to the appropriate directory within the data/ref subdirectory as part of the setup.

4. One important note is that after the STAR alignment step, I moved the unmapped reads to a separate directory (data/unmapped) for ease of access and future reference. This step is mentioned in the protocol_README.md file.

5. Contigs < 300 bp were filtered after BLAST analysis in R and individual excel files for each virus identified were created using an R script.

6. The R directory was uploaded separately and it contains seperate excel files for each virus identified in the analysis. The R script used to generate these files is included in the scripts directory as virus_script.R

7. All files in the FP directory are relevant to the final project submission.




