# DAP-seq data processing pipeline

## Pipeline structure

The DAP-seq pipeline consists of three major steps: sequence reads 
preprocessing, mapping of the reads to reference genome and peak calling.
Additional postprocessing steps include motif search and mapping of 
peaks to coding genes.

Separate shell scripts run each step of the pipeline. You can inspect 
results of each step before running the next script. 

## Requirements

The pipeline runs on OS Linux (developed and tested on Ubuntu 16.04). 

Some external programs should be installed before running the pipeline:

* SolexaQA++ (http://solexaqa.sourceforge.net/)
* FOCUS (https://sourceforge.net/projects/metagenomefocus/)
* bowtie, v.1.1.2 (http://bowtie-bio.sourceforge.net/index.shtml)
* samtools (http://www.htslib.org/)
* MACS2 (https://github.com/taoliu/MACS)
* MEME Suite (http://meme-suite.org/)

bowtie, samtools and MACS2 can be installed from Ubuntu package manager. 
Other programs should be linked to /usr/bin or other directory,
 where binaries are stored normally. Pipeline depends on "SolexaQA++", 
 "focus.py" and "meme" commands.

## Reference data

Out of the box, three reference genomes are available: Pseudomonas stutzeri
RCH2, Pseudomonas fluorescens FW300-N2E2 and P. putida KT2440. 
To build Bowtie indices for these reference genomes, execute script 
build_bowtie_indices.sh before running the pipeline.

## Input files

The pipeline works with single-end sequence reads in gzipped FASTQ file. 
All FASTQ files must be in the same directory. 
All FASTQ files must have names XX_DAP_NN_R1_001.fastq.gz, where XX is
genome label (any letters or numbers), and NN is the sample number 
(any number of any digits). For example: RCH2_DAP_012_R1_001.fastq.gz.
One additional file (called "directory file") with descrption of samples 
required for configuration the pipeline. The directory file contains six
tab-separated fields for each sample in the following order:

* Sample label (FASTQ file name without ending "_R1_001.fastq.gz")
* Reference genome name (as defined in perl/libs.tsv).
* Protein ID
* Replicate number (1, 2, 3...)
* Treatment label (like "Phosphorylated", "N/A" etc.)
* Label of control sample (one of labels from column 1. It will be used for peak calling. Usually, all samples have the same control label)

Directory file may have comment lines starting with "#".

Example:

|#Sample|Genome|Protein|Replicate|Treatment|Control|
|-------|------|-------|---------|---------|-------|
|N2E2_01|P. fluorescens N2E2|Pf6N2E2_1296|1|N/A|N2E2_99|


## How to generate pipeline

To create the shell scripts running three steps of the pipeline, run
 configuration Perl scripts configure_dapseq_qc_mapping.pl and configure_dapseq_peak_calling.pl. 
If you want to run motif search, run additional configuration script 
configure_meme.pl

## Configuration of read preprocessing and mapping

Run
```
cd perl
perl configure_dapseq_qc_mapping.pl <full path to directory file> <path to directory with FASTQ files> <working directory> <output directory for BAM files> <output directory for FOCUS files>
```

This script creates two shell scripts **in the working directory**: run_qc.sh 
and run_mapping.sh

The run_qc.sh script runs SolexaQA++ for read trimming and removal of 
reads shorter than 30 bp. 

For the reads that pass read trimming, FOCUS generates phylogenetic profile. 
If you have many unmapped reads, inspect FOCUS output 
for possible cross-contamination with DNA from unrelated organisms.

Output sequence files are gzipped.

The run_mapping.sh script runs bowtie (v.1) for read mapping to the reference 
genome. 
