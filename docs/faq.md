---
layout: default
title: FAQ
---
##Frequently Asked Questions

### Can you make the re-assembly of the Wrighton data available?

The re-assembly of the [Wrighton et al., 2012](http://www.sciencemag.org/content/337/6102/1661) data can be found in the [Additional Data](additional.html) section.

### I can not match your sequences on NCBI with the data in your script?

The sequences on NCBI are contigs. They have a header structure that can be used to relate them to the parrent scaffold. E.g. HPminus3.2 is the second contig in scaffold number 3. **Scaffold number 3** is the ID used in the R data to associate e.g. coverage to the scaffold. 

Alternatively you can download the scaffolds [here](https://www.dropbox.com/s/ig6gie43tqfziif/assembly.rar).

### I do not have access to CLC - how do I map reads to the assembly and generate the .csv file used in R?

There are plenty of free great short read mapping software available. We usually use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) if we need an alternative to CLC. To map reads to an assembly you first need to generate an indexed version of the de novo assembly (`assembly.fa`) using `bowtie2-build`.

{% highlight text%}
bowtie2-build assembly.fa assembly
{% endhighlight %}

You can now map the reads to the assembly and store the alignment in SAM format using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

{% highlight text%}
bowtie2 -x assembly -U reads.fastq -S alignment.sam
{% endhighlight %}

From the alignment (`alignment.sam`) there are many ways to extract the coverage of each scaffold. One option is to use the `depth` function in [samtools](http://samtools.sourceforge.net/).

Using [samtools]((http://samtools.sourceforge.net/) first convert the `.sam` file to `.bam` format and then sort it. Now the `depth` function is used to generate a position based coverage profile of each scaffold (`depth.txt`). 

{% highlight text%}
samtools view -bS alignment.sam > alignment.bam
samtools sort alignment.bam alignment.sorted
samtools depth alignment.sorted.bam > depth.txt
{% endhighlight %}

To generate single average coverage vaule for each scaffold you can use the small perl script `calc.coverage.in.bam.depth.pl` in the folder `multi-metagenome/misc.scripts`. It generates `coverage.length.csv.` which is equivalent to the `.csv` file obtained using CLC.

{% highlight text%}
perl calc.coverage.in.bam.depth.pl -i depth.txt -o coverage.length.csv
{% endhighlight %}

Alternative you can load `depth.txt` directly in R and make the calculations in there. However, the size of the `depth.txt` file might be quite large as it stores a line for each base in the assembly. Anyway if you want to handle the data in R you could do something like this:

{% highlight r%}
cov <- read.delim("depth.txt", header = F)
colnames(cov) <- c("scaffold", "position", "coverage")
cov.scaffold <- as.data.frame(tapply(cov$coverage, cov$scaffold, mean))
length.scaffold <- as.data.frame(tapply(cov$position, cov$scaffold, max))
{% endhighlight %}

### I do not have access to CLC - can you reccomend other metagenome assemblers?

There are several free de novo assembly programs that are able to handle metagenome data. See e.g. [SOAPdenovo](http://soap.genomics.org.cn/soapdenovo.html), [metAMOS](http://cbcb.umd.edu/software/metAMOS), [Ray Meta](http://denovoassembler.sourceforge.net/) or [IDBA-UD](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/index.html).

Alternatively the program [khmer](https://khmer.readthedocs.org/en/latest/) can be used as a preprocessing step to enable de novo assembly of very large datasets.



