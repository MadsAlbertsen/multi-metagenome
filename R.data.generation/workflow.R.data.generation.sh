#!/bin/bash

##### Description
# This small shell script wraps the data generation for R in Albertsen et. al., 2013

##### Needed input files
# Fasta file with all assembled scaffolds (keep the naming as >1, >2 etc): assembly.fa

##### Needed software
# Prodigal
# HMMER 3.0
# BLAST
# MEGAN
# Perl scripts from : git clone git://github.com/MadsAlbertsen/multi-metagenome.git        

clear
echo "---Metagenomics workflow script v.1.0---"

echo ""
echo "Calculating tetranucleotide frequency"
perl multi-metagenome/R.data.generation/calc.kmerfreq.pl -i assembly.fa -o assembly.kmer.tab

echo ""
echo "Calculating gc content"
perl multi-metagenome/R.data.generation/calc.gc.pl -i assembly.fa -o assembly.gc.tab

echo ""
echo "Finding essential genes - Predicting proteins (Prodigal)"
prodigal -a temp.orfs.faa -i assembly.fa -m -o temp.txt -p meta -q
cut -f1 -d " " temp.orfs.faa > assembly.orfs.faa

echo ""
echo "Finding essential genes - running HMM search"
hmmsearch --tblout assembly.hmm.orfs.txt --cut_tc --notextw multi-metagenome/R.data.generation/essential.hmm assembly.orfs.faa > hmm.temp.txt
tail -n+4  assembly.hmm.orfs.txt | sed 's/ * / /g' | cut -f1,4 -d " " | sed 's/_/ /' > assembly.orfs.hmm.id.txt
grep -v "#" assembly.hmm.orfs.txt | cut -f1 -d " " > list.of.positive.orfs.txt
perl multi-metagenome/R.data.generation/extract.using.header.list.pl -l list.of.positive.orfs.txt -s assembly.orfs.faa -o assembly.orfs.hmm.faa

echo ""
echo "Finding essential genes - Blasting positive hits"
blastp -query assembly.orfs.hmm.faa -db refseq_protein -evalue 1e-5 -num_threads 60 -max_target_seqs 5 -outfmt 5 -out assembly.orfs.hmm.blast.xml

echo ""
echo "Finding essential genes - Extracting consensus taxonomic assignment"
MEGAN +g -x "import blastfile= assembly.orfs.hmm.blast.xml meganfile=temp.rma;recompute toppercent=5;recompute minsupport=1;update;collapse rank=Species;update;select nodes=all;export what=CSV format=readname_taxonpath separator=tab file=assembly.orfs.hmm.blast.tax.txt;update;close"
perl multi-metagenome/R.data.generation/hmm.majority.vote.pl -i assembly.orfs.hmm.blast.tax.txt -o assembly.tax.consensus.txt -n
sed 's/\t/;/' assembly.orfs.hmm.blast.tax.txt | cut -f1,5 -d ";" | sed 's/;/\t/' | sed 's/_/\t/'  > assembly.orfs.hmm.blast.tax.tab

#echo ""
#echo "Removing temp files"
#rm hmm.temp.txt
#rm list.of.positive.orfs.txt
#rm assembly.orfs.hmm.blast.tax.txt
#rm temp.orfs.faa
#rm temp.txt
#rm temp.rma
#rm assembly.orfs.hmm.blast.xml
#rm assembly.orfs.hmm.faa
