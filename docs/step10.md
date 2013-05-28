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
The only file needed is a [SAM file](http://samtools.sourceforge.net/) of the read mappings to the assembled scaffolds, which is produced by most read mapping tools. The SAM file contains all information of the allignment of each specific read to the assembled scaffolds. Note that the script is very simple in it's approach to identifying PE reads. It assumes a read naming structure of `read1_1` and `read1_2` which is the default structure of the raw Illumina reads. A full readname of a Illumina PE read pair can be seen below, note that the readname prior to the underscore is identical between the two reads.

{% highlight text%}
HWI-ST1040:48:c06lnacxx:2:1101:7493:111092_1:N:0:CAGATC
HWI-ST1040:48:c06lnacxx:2:1101:7493:111092_2:N:0:CAGATC
{% endhighlight %}

###Parameters 
The script `cytoscapeviz.pl` searches through the SAM files and identifies reads that map to different scaffolds. The output is a file that states which scaffolds that are linked together. There are a few options that can be used to control the output. 

`f` controls the minimum number of connections between scaffolds needed to call a link. In the example below at least 2 PE links are needed before the scaffolds are reported to be linked. 

`e` is used to only select links that map to the ends of the scaffolds, `e` = 500 means that only PE links mapping within 500 bp of each end is used to call links between scaffolds. 

`m` is used to set a minimum length of scaffolds that are reported to be circular, i.e. where a PE read pair map to each end of the same scaffold. 

`a` is the average read length of the reads and is used to calculate a coverage for each scaffold (It could also be generated directly from the SAM file, but this simple speeds up the script a little). 

`c` is a flag to get an additional condensed output. By default `cytoscapeviz.pl` generates connections for each end of all scaffolds. E.g. `scaffold 1` will be represented by the nodes `scaffold 1 start` and `scaffold 1 end`, however this can be a little messy to look at initially. Instead `c` is a flag stating that an additional condensed output should be made, where all scaffolds are represented by a single node (e.g. `scaffold 1`) instead of a start and end node. Always start by looking at the condensed output.

{% highlight text%}
perl \multi-metagenome\cytoscapeviz\cytpscapeviz.pl -i mapping.sam -f 2 -e 500 -m 3000 -a 125 -c
{% endhighlight %}

###Output



##Visualising PE connections using Cytoscape
Describe how to use cytoscape and show a few use cases
