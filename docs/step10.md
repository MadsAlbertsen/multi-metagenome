---
layout: default
title: PE tracking
---
#Under Construction
##Tracking paired-end reads
As some genes are present in multiple copies (e.g. 16S rRNA genes or transposases) they will not be included in the initial coverage-defined subset (see [Binning](step9.html)). Instead paired-end read information is used to associate multiple copy genes with the appropriate genome bin.

##Background
Paired-end (PE) reads are often generated when using the Illumina sequencing platform. DNA is fragmented into smaler pieces (200-700 bp) and then sequenced from each end. As the DNA fragments are larger than the length that can be sequenced, there will be unknown sequence between the two reads. However, the two reads originate from the same piece of DNA and can therefore be considered linked to one another.

![PE definition](figure/PEdef.png)

Paired-end information is utilized in the process of scaffolding, where assembled contigs are merged into scaffolds based on PE links. Here `contig A` and `contig B` can be merged together into `scaffold AB` as they are connected by a PE read.

![PE connection](figure/PEcon.png)

Now consider the case where we have a repeat in the genome (I.e. a piece of DNA in multiple copies that can't be spanned by PE reads). Here we can not make a scaffold as we do not know if it should be `scaffold CE` or `scaffold CF` and `scaffold DE` or `scaffold DF`. However, we do know that the repeat is associated with contigs `E`, `C`, `D` and `F` and can therefore be said to originate from the same genome.

![PE repeat](figure/PErepeat.png)

It's not just repeats that can be identified in this way. Short contigs can by chance (or extreme GC content) have a coverage that is significant different from the rest of the genome and therefore not included in the initial coverage defined bin. These can normally be picked up in through tracking of PE reads.

##Generating PE connections
Describe cytoscapeviz.pl

##Visualising PE connections
Describe how to use cytoscape and show a few use cases
