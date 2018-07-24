# How to add new refrence genome

This page describes how to add new genome to the reference library. Currently, 
DAP-seq pipeline works only with organisms from [Fitness Browser](http://fit.genomics.lbl.gov/). 
See list of available orgenisms [here](http://fit.genomics.lbl.gov/cgi-bin/orgAll.cgi).

## What you need

* Genome sequences in FASTA format (one file with all replicons)
* Genome sequences in GenBank format (separate files for each replicons)
* List of predicted operons (tab-separated). Each line must contain operon ID, start position, end position, strand, number of genes, list of genes (use [data/rch2/operons.tsv] as an example).
* List of genes with functional annotations (tab-separated). Each line must contain gene ID, locus tag, gene function, genome ID, genome name (use [data/rch2/img.tsv] as an example).
* List of genes from [Fitness Browser](http://fit.genomics.lbl.gov/cgi-bin/orgAll.cgi)

## What to do

1. Create new directory in the data subfolder and copy all your files in it.
2. Create bowtie index with your FASTA file as input (consult manual 
   for bowtie-build command for details).
3. Create file "replicons.tsv" in the same folder with list of all GenBank files and sequence lengths (use [data/rch2/replicons.tsv] as an example).
3. Edit data/libs.tsv file. Add a new line with eight fields (tab-separated):
      1. genome name
      2. relative path to bowtie index
      3. effective genome size (to be used by MACS2)
      4. relative path to file with gene functions from [Fitness Browser](http://fit.genomics.lbl.gov/cgi-bin/orgAll.cgi)
      5. name of your genome in [Fitness Browser](http://fit.genomics.lbl.gov/cgi-bin/orgAll.cgi) database
      6. relative path to file with gene functions
      7. relative path to file with predicted operons
      8. relative path to file with list of replicons
4. Edit data/genome_list.txt file. Add a new line with four fields (tab-separated)
   for each of your GenBank files:
      1. NCBI accession (with version)
      2. Genome name (as in libs.tsv)
      3. relative path to GenBank file.
      4. "chr" for chromosomes or "plasmid" for plasmids

## What if you organism is missing from Fitness Browser

Create empty file for list of genes from Fitness Browser. Leave blank the field for 
name of organism in Fitness Browser. You will still be able to run preprocessing,
mapping and peak calling, but not the report generation.
