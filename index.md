---
layout: default
title: index
---
This project contains scripts and tutorials on how to assemble individual microbial genomes from metagenomes, as described in:

**Genome sequences of rare, uncultured bacteria obtained by differential coverage binning of multiple metagenomes**

[Mads Albertsen](http://personprofil.aau.dk/120257), [Philip Hugenholtz](http://ecogenomic.org/users/phil-hugenholtz), [Adam Skarshewski](http://ecogenomic.org/users/adam-skarshewski), [Gene W. Tyson](http://www.ecogenomic.org/users/gene-tyson), [Kåre L. Nielsen](http://personprofil.aau.dk/103057) and [Per .H. Nielsen](http://personprofil.aau.dk/105842)

Nature Biotechnology 2013, doi: [10.1038/nbt.2579](http://www.nature.com/nbt/journal/vaop/ncurrent/abs/nbt.2579.html)

##What is differential coverage binning?
Differential coverage binning refers to the use of the abundance of the bacteria in the samples as the primary method of extracting genomes from metagenomes. This is compared to many other binning methods where sequence composition (e.g. tetranucleotide patterns) is used for binning, which is hampered by local sequence deviations within genomes. 

The great advantage of using abundance is that it allows us to use abundance estimates from multiple samples, thereby increasing the binning resolution greatly. Given that sequencing prices continiues to drop it is already now cheaper to generate data than to analyse it.

##Step-by-step guide
The [guide](docs/overview.html) covers all aspects, from which samples that can be used, to binning, finishing and validation of the extracted genomes. [Rstudio](http://www.rstudio.com/) (a powerfull IDE to [R](http://www.r-project.org/)) is used as the main tool for data handling as it allows integration of all relevant data, which is key for dealing with large and complex datasets. The guide to binning in R is also available in R markdown format [here](https://github.com/MadsAlbertsen/multi-metagenome/tree/master/R.markdown.guide), which allows direct recreation of the binning steps in R from the original processed data used in the paper.

In the [overview](docs/overview.html) section a short description of the workflow is given, along with a detailed workflow figure that summarise the different steps in the process of assembling individual genomes from metagenomes.

In addition to the online guide a PDF version of the original published guide is available [here](https://github.com/MadsAlbertsen/multi-metagenome).

##Data availability
If you want to reanalyze the data used in this study, the raw fastq reads can be obtained from NCBI SRA: [SRX206471 (HP+)](http://www.ncbi.nlm.nih.gov/sra/SRX206471?report=full) and [SRX247688 (HP-)](http://www.ncbi.nlm.nih.gov/sra/SRX247688?report=full). The assembled contigs can be obtained from NCBI GenBank under accesion number [APMI01000000](http://www.ncbi.nlm.nih.gov/nuccore/494587257). 

In addition, all processed data ready for binning in R can be obtained [here](https://dl.dropbox.com/s/989dix16ugyuvrq/Albertsen2013.data.tar.gz).

