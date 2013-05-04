# Metagenome assembly guide

This project contains scripts and tutorials on how to assemble individual microbial genomes from metagenomes, as described in Albertsen *et al*., 2013 (add link):

**Recovery of genomes from rare, uncultured bacteria by differential coverage binning of multiple deep metagenomes**

##  Introduction
Some text

## Data availability
If you want to reanalyze the data used in this study the raw fastq reads can be obtained from the NCBI short read archive [SRX206471 (HP+)](http://www.ncbi.nlm.nih.gov/sra/SRX206471?report=full) and [SRX247688 (HP-)](http://www.ncbi.nlm.nih.gov/sra/SRX247688?report=full). The assembled contigs can be obtained from NCBI GenBank under accesion number APMI01000000 (note: link not updated by NCBI yet). 

In addition all processed data used in R can be obtained from HERE (add link).

## Step-by-step guide
A detailed step-by-step guide is available in [this folder](https://github.com/MadsAlbertsen/multi-metagenome), which explains all steps from sampling to extraction of individual genomes.

A small shell script can be found in [R.data.generation] that wraps the individual scripts needed for generation of all data needed for subsequent genome extraction through R. The only input needed is a fasta file of the assembled scaffolds.

The R part of the step-by-step guide is also available as a seperate tutorial in R markdown format, which generates the main figures used in the article. Check it out HERE (add link).

05/2013 - Mads Albertsen - contact: ma at bio dot aau dot dk
