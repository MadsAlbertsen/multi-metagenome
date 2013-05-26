---
layout: default
title: Assembly
---
##Sequencing
We normally sequence `2x150 bp reads` on the Illumina HiSeq2000 platfrom, using the nextera paired-end protocol (PE, 200-700 bp insert size). While PE libraries are enough for the assembly itself and enables recovery of the complete genome, mate-pair libraries (MP, >2 kbp insert size) are extremely valuable in order to obtain the genomes in a single chromosome. The vast majority of repeats in bacterial genomes are due to transpoases (1-1.5 kbp) which can be spanned by MP libraries.

###How much sequencing is needed?
It is straight forward to get a ball-park estimation of the sequencing depth needed to recovere your favorite organism if you know the approximate abundance of it. As an example we assume an average genome size of 3 Mbp of all bacteria in the sample and the abundance of the target organism to 1%. In order to obtain a good draft genome at least 50X coverage is needed using 150 bp PE Illumina data. Hence, the needed sequencing depth would be:
{% highlight text%}
Genomesize * Coverage * Abundance = Needed sequencing depth
{% endhighlight %}
{% highlight r%}
3 Mbp * 50 X * 1 % = 15 Gbp 
{% endhighlight %}
15 Gbp is currently 1/2 lane of HiSeq2000 data or 2 MiSeq runs. Note that the needed sequencing depth can be split between the different samples.

##Assembly
We use [CLC genomics workbench](http://www.clcbio.com/) for all trimming, assembly and read mapping out of convinience, but all parts can be substituted with open source tools. 

###Trimming
Prior to assembly the reads needs to be trimmed in order to remove adapters and bad quality sequence. We normally apply the following criteria:

{% highlight text%}
Remove Illumina nextera adapters if found
A minimum phred score of 20
Allowing no ambigious nucleotids (N's)
A minimum read length of 50 bp
{% endhighlight %}

###De novo assembly
De novo assembly of metagenomes is a field in constant development and numerous strategies exists. We are currently using CLC's de novo assembly implementation as it is able to handle a wide range of coverage abundances, is memory friendly and fast. We have sucessfully assembled metagenome datasets >300 Gbp in less than a day on server with 40 processers and 256 Gbp of RAM. We use standard settings except a `kmer of 63`.

If CLC is not an option (It's quite expensive) we recommend to take a look at [khmer](https://khmer.readthedocs.org/en/latest/) which also tries to enable assembly of large metagenome datasets.

However, the main problem in metagenome assembly is the presence of multiple closely realted strains, which acts like poison to the assemblers. Although several groups are trying to tackle this, there is to our knowledge no good bioinformatic solution. 

[Next: Data generation](step5.html)
