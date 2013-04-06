#!/usr/bin/env perl
###############################################################################
#
#    circosviz.pl version 1.0
#    
#    Indicates connections between contigs/scaffolds and produces circos 
#    files for visualization 
#
#    Copyright (C) 2012 Mads Albertsen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;
use POSIX;            #to use the floor command

#core Perl modules
use Getopt::Long;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params

my $global_options = checkParams();

my $samfile;
my $fastafile;
my $enddist;
my $minlength;
my $avgreadlength;
my $headersplit;
my $circoscbin;
my $minpesplit;
my $frc;

$samfile = &overrideDefault("data.sam",'samfile');
$fastafile = &overrideDefault("data.fasta",'fastafile');
$enddist = &overrideDefault("500",'enddist');
$minlength = &overrideDefault("3000",'minlength');
$avgreadlength = &overrideDefault("126",'avgreadlength');
$headersplit = &overrideDefault("_",'headersplit');
$circoscbin = &overrideDefault("10000",'circoscbin');
$minpesplit = &overrideDefault("1000",'minpesplit');
$frc = &overrideDefault("N",'frc');

my $line;
my $old;
my $new;
my $count = 0;
my $connectionorder;
my $readcount = 0;
my $printreadcount = 0;
my $dummy = 0;
my $seq2;
my $contigcolor;
my $bandstart;
my $bandend;
my $firstn = 0;
my $contiglength;
my $bandcount;
my $countcends = 0;
my $countcdiff = 0;
my $countcsamewl = 0;
my $numcontigs = 0;
my $totallength = 0;
my @splitline;
my @splitline1;
my @splitline2;
my @contigname;
my @contiglength;
my @seqarray;
my %contigs;
my %contigcov;
my %reads;
my %breads;
my %connections;
my %contigswhit;
my %circoscov;
my %circoscovE;
my %circoscovD;
my %circoscovS;
my %counts;
my %karyotypeh;

######################################################################
# CODE HERE
######################################################################

open(INsam, $samfile) or die("Cannot open $samfile\n");
open(INfasta, $fastafile) or die("Cannot open $fastafile\n");
open(OUTcends, ">circos.ends.txt") or die("Cannot create circos.ends.txt\n");   
open(OUTcdiff, ">circos.dcontigs.txt") or die("Cannot create circos.dcontigs.txt\n");   
open(OUTcsamewl, ">circos.scontigs.wl.txt") or die("Cannot create circos.scontigs.wl.txt\n");   
open(OUTckary, ">circos.karyotype.txt") or die("Cannot create circos.karyotype.txt\n");   
open(OUTccov, ">circos.coverage.txt") or die("Cannot create circos.coverage.txt\n");   
open(OUTccovE, ">circos.count.ends.txt") or die("Cannot create circos.coverageE.txt\n");   
open(OUTccovD, ">circos.count.dcontigs.txt") or die("Cannot create circos.coverageD.txt\n");   
open(OUTccovS, ">circos.count.scontigs.wl.txt") or die("Cannot create circos.coverageS.txt\n");   
open(OUTcgc, ">circos.gc.txt") or die("Cannot create circos.gc.txt\n");   
open(OUTcrules, ">circos.rules.txt") or die("Cannot create circos.rules.txt\n");   

if (my $is_even = $avgreadlength % 2 == 0){                                                                               #To be able to get a better estimate of where the read maps in the start and end / otherwise it would just have used the start position
	$avgreadlength = $avgreadlength/2;
}
else{
	$avgreadlength = ($avgreadlength+1)/2;
}

print "Generating overview of pe connections in the SAM file.\n";

while ( $line = <INsam> ) {
	chomp $line;                                                         												     
	$count++;
	@splitline = split(/\t/, $line); 	
	if ($line =~ m/\@SQ/) {                               												                  #if we are in the contig header area then retrive all contigs/scaffolds and store the name and length in the hash: contig		
			@contigname = split(/:/, $splitline[1]);                     												  #Retrive the contig name
			@contiglength = split(/:/, $splitline[2]);																	  #Retrive the contig length
			$contigs{$contigname[1]} = $contiglength[1];                   												  #Make a hash with key = "contig name" and value = "contig length"			
			$totallength = $totallength + $contiglength[1];                   											
			$contigcov{$contigname[1]} = 0;
			$numcontigs++;
		}	
	else {
		if ($line !~ m/(\@PG|\@HD|\@SQ)/) { 
			@splitline = split(/\t/, $line);
			my $covbin = floor($splitline[3]/$circoscbin);
			$contigcov{$splitline[2]}++;	
			my $tempcovbin = $splitline[2].":".$covbin;
			$circoscov{$tempcovbin}++;	
			@splitline1 = split(/$headersplit/, $splitline[0]);															  #The read header - assumes the new Illumina format! e.g. "name space 1" or "name space 2"		
			if ($splitline[1] != 19 and $splitline[1] != 35){															  #SAM Flags that indicate that these PE reads are maping as they are supposed to.. hence they are not interesting..																									#Between diferent contigs? Between the same contig - e.g check if plasmid.
				if ((($splitline[3]+$avgreadlength) <= $enddist or ($splitline[3]+$avgreadlength) >= ($contigs{$splitline[2]}-$enddist))  and !exists($breads{$splitline1[0]})) {            #The read is required to hit within a certain distance from the contigs ends.. The middle postition of the read is used
					if (exists($reads{$splitline1[0]})){                      											  #If one of the PE reads has already been seen then add the hit to the match hash														
						@splitline2 = split(/\t/,$reads{$splitline1[0]});	
						delete $reads{$splitline1[0]};
						if ($splitline[2] ne $splitline2[0]){                                                             #Good connection diff contigs
							my $tempend = $splitline[3]+1;
							print OUTcends "$count $splitline[2] $splitline[3] $tempend\n";
							my $tempcovbin2 = $tempcovbin;
							if ($covbin != 0){
									$tempcovbin2 = $splitline[2].":1";
								}
							$circoscovE{$tempcovbin2}++;
							$tempend = $splitline2[1]+1;	
							print OUTcends "$count $splitline2[0] $splitline2[1] $tempend\n";	
							$tempcovbin2 = $splitline2[0].":0";
							if (floor($splitline2[1]/$circoscbin) != 0){
								$tempcovbin2 = $splitline2[0].":1";
							}
							$circoscovE{$tempcovbin2}++;								
							$countcends++;
							
						}
						else{			                                                                                  
							if ($contigs{$splitline2[0]} >= $minlength and abs($splitline[3]-$splitline2[1])>= $minpesplit){   #Good circular connection								
								my $tempend = $splitline[3]+1;
								print OUTcends "$count $splitline[2] $splitline[3] $tempend\n";
								my $tempcovbin2 = $tempcovbin;
								if ($covbin != 0){
									$tempcovbin2 = $splitline[2].":1";
								}
								$circoscovE{$tempcovbin2}++;
								$tempend = $splitline2[1]+1;
								print OUTcends "$count $splitline2[0] $splitline2[1] $tempend\n";
								$tempcovbin2 = $splitline2[0].":0";
								if (floor($splitline2[1]/$circoscbin) != 0){
									$tempcovbin2 = $splitline2[0].":1";
								}
								$circoscovE{$tempcovbin2}++;									
								$countcends++;
							}
						}
					}
					else{								
						$reads{$splitline1[0]} = "$splitline[2]\t$splitline[3]";										   #If the other PE read has not been seen then create the first instance of the pair
					}
				}
				else{
					#Not end connection but still different length
					@splitline1 = split(/$headersplit/, $splitline[0]);														#The read header - assumes the new Illumina format! e.g. "name space 1" or "name space 2"
					if (exists($reads{$splitline1[0]}) or exists($breads{$splitline1[0]})){                      			 #If one of the PE reads has already been seen then add the hit to the match hash														
						if (exists($reads{$splitline1[0]})){						
							@splitline2 = split(/\t/,$reads{$splitline1[0]});												#get the contig name from the old read						
							delete $reads{$splitline1[0]};
						}
						else{
							@splitline2 = split(/\t/,$breads{$splitline1[0]});
							delete $breads{$splitline1[0]};
						}
						if ($splitline[2] ne $splitline2[0]){                                                               #bad connection diff contigs
							my $tempend = $splitline[3]+1;
							print OUTcdiff "$count $splitline[2] $splitline[3] $tempend\n";
							$circoscovD{$tempcovbin}++;
							$tempend = $splitline2[1]+1;
							print OUTcdiff "$count $splitline2[0] $splitline2[1] $tempend\n";
							my $tempcovbin2 = $splitline2[0].":".floor($splitline2[1]/$circoscbin);
							$circoscovD{$tempcovbin2}++;	
							$countcdiff++;			
						}
						else{			                                                                                  
							if ($contigs{$splitline2[0]} >= $minlength and abs($splitline[3]-$splitline2[1])>= $minpesplit){   #bad same connection								
								my $tempend = $splitline[3]+1;
								print OUTcsamewl "$count $splitline[2] $splitline[3] $tempend\n";
								$tempend = $splitline2[1]+1;
								$circoscovS{$tempcovbin}++;
								print OUTcsamewl "$count $splitline2[0] $splitline2[1] $tempend\n";
								my $tempcovbin2 = $splitline2[0].":".floor($splitline2[1]/$circoscbin);
								$circoscovS{$tempcovbin2}++;
								$countcsamewl++;
								
							}
						}
					}
					else{								
						$breads{$splitline1[0]} = "$splitline[2]\t$splitline[3]";										    #If the other PE read has not been seen then create the first instance of the pair -
					}
				}
			}			
			else{																											#here there is room for adding some general contig stats - e.g. the proportion of PE to SE reads? or the contig coverage stats																	
			}
			$readcount++;                                                                                                   #keep track of the number of reads look though
			$printreadcount++;
			if ($printreadcount == 1000000) {
				$printreadcount = 0;
				print "$readcount reads looked through\n";
			}		
		}	
	}
}
close INsam;

print "Generating coverage file.\n";


my @tempcov = keys %circoscov;
@tempcov = sort @tempcov;
foreach my $cov (@tempcov){
	my @splitline = split(/:/,$cov);	
	my $start = $splitline[1]*$circoscbin;
	my $end = ($splitline[1]+1)*$circoscbin;
	if ($end > $contigs{$splitline[0]}){
		$end = $contigs{$splitline[0]};
	}
	my $tempcov = $circoscov{$cov}/($end-$start)*$avgreadlength*2;                                                          #*2 since I divided it with 2 in the beginning.. 
	print OUTccov "$splitline[0] $start $end $tempcov\n"; 
}

@tempcov = keys %circoscovE;
@tempcov = sort @tempcov;
foreach my $cov (@tempcov){
	my @splitline = split(/:/,$cov);	
	my $start = $splitline[1]*$circoscbin;
	my $end = $enddist;
	if ($splitline[1] == 1){
		$start = $contigs{$splitline[0]} - $enddist;
		$end = $contigs{$splitline[0]};		
	}		
	print OUTccovE "$splitline[0] $start $end $circoscovE{$cov}\n"; 
}

@tempcov = keys %circoscovD;
@tempcov = sort @tempcov;
foreach my $cov (@tempcov){
	my @splitline = split(/:/,$cov);	
	my $start = $splitline[1]*$circoscbin;
	my $end = ($splitline[1]+1)*$circoscbin;
	if ($end > $contigs{$splitline[0]}){
		$end = $contigs{$splitline[0]};
	}
	print OUTccovD "$splitline[0] $start $end $circoscovD{$cov}\n"; 
}

@tempcov = keys %circoscovS;
@tempcov = sort @tempcov;
foreach my $cov (@tempcov){
	my @splitline = split(/:/,$cov);	
	my $start = $splitline[1]*$circoscbin;
	my $end = ($splitline[1]+1)*$circoscbin;
	if ($end > $contigs{$splitline[0]}){
		$end = $contigs{$splitline[0]};
	}
	print OUTccovS "$splitline[0] $start $end $circoscovS{$cov}\n"; 
}



print "Generating karyotype file based on fasta file.\n";

while (my $line = <INfasta>)  {
	if ($line =~ m/>/) {
		chomp $line;
		$line =~ s/\>//g; 	
		if ($dummy == 1){
			push (@seqarray, $seq2);
			@splitline = split(/\t/,$seq2);
			$contiglength = length($splitline[1]);
			$karyotypeh{$splitline[0]}="chr - $splitline[0] $splitline[0] 0 $contiglength set1-7-qual-";
		}
		$seq2 = "$line\t";
		$dummy =1;
	}
	else {
		chomp $line;
		$seq2 = $seq2.$line;
		}
}
close INfasta;
push (@seqarray, "$seq2");	          #to catch the last sequence
@splitline = split(/\t/,$seq2);
$contiglength = length($splitline[1]);
$karyotypeh{$splitline[0]}="chr - $splitline[0] $splitline[0] 0 $contiglength set1-7-qual-";                                #Change the color scheme here if needed
foreach my $key (reverse sort {$contigs {$a} <=> $contigs {$b}} keys %contigs){
	$contigcolor++;
	if ($contigcolor == 8){																			    			        #Change the color scheme here if needed
		$contigcolor = 1;
	}
	print OUTckary "$karyotypeh{$key}$contigcolor\n";
	print OUTcrules "<rule>\n";
	print OUTcrules "condition = _CHR1_ eq ",'"',"$key",'"',"\n";
	print OUTcrules "color = set1-7-qual-$contigcolor",'_a5',"\n";                                                          #Change the color scheme here if needed
	print OUTcrules "</rule>\n";
	
}
print OUTckary "\n";                                                                                                        #Just to seperate for the band part of the karyotype file



print "Generating gc file based on fasta file.\n";

foreach my $sequence (@seqarray){
	$counts{G} = 0;
	$counts{C} = 0;
	$counts{A} = 0;
	$counts{T} = 0;
	$firstn = 0;
	my @tseq = split("\t", $sequence);
	my @seq = split("", $tseq[1]);
	my $bincount = 0;
	my $lengthcount = 0;
	my $lastend = 0;
	foreach my $nucleotide (@seq) {
		$counts{$nucleotide}++;
		$bincount++;
		$lengthcount++;		
		if ($bincount == $circoscbin){
			my $gc = ($counts{G}+$counts{C})/($counts{G}+$counts{C}+$counts{A}+$counts{T})*100;
			print OUTcgc "$tseq[0] $lastend $lengthcount $gc\n";
			$lastend = $lengthcount;
			$counts{G} = 0;
			$counts{C} = 0;
			$counts{A} = 0;
			$counts{T} = 0;
			$bincount = 0;
		}
		if ($nucleotide eq "N"){
			if ($firstn == 0){
				$firstn = 1;
				$bandstart = $lengthcount;
			}				
		}
		else{
			if ($firstn == 1){
				$firstn = 0;
				$bandcount++;
				print OUTckary "band $tseq[0] band$bandcount band$bandcount $bandstart $lengthcount black\n";
			}
		}
	}
	my $gc = ($counts{G}+$counts{C})/($counts{G}+$counts{C}+$counts{A}+$counts{T})*100;
	print OUTcgc "$tseq[0] $lastend $lengthcount $gc\n";	
}
 

if ($frc ne "N"){
	open(INfrc, $frc) or die("Cannot open $frc\n");
	open(OUTfrc, ">frc.tracks.txt") or die("Cannot create frc.tracks.txt\n");   
	while ( my $line = <INfrc> ) {
		chomp $line;   	
		next if ($line =~ m/#/);
		my @splitline = split(/\t/,$line);
		if ($line =~ m/STRECH_PE/) {$splitline[2] = 1};
		if ($line =~ m/COMPR_PE/) {$splitline[2] = 2};
		if ($line =~ m/LOW_NORM_COV_PE/) {$splitline[2] = 3};
		if ($line =~ m/LOW_COV_PE/) {$splitline[2] = 4};
		if ($line =~ m/HIGH_COV_PE/) {$splitline[2] = 5};
		if ($line =~ m/HIGH_NORM_COV_PE/) {$splitline[2] = 6};
		if ($line =~ m/HIGH_OUTIE_PE/) {$splitline[2] = 7};
		if ($line =~ m/HIGH_SINGLE_PE/) {$splitline[2] = 8};
		if ($line =~ m/HIGH_SPAN_PE/) {$splitline[2] = 9};
		if ($line =~ m/COMPR_MP/) {$splitline[2] = 10};
		if ($line =~ m/HIGH_OUTIE_MP/) {$splitline[2] = 11};
		if ($line =~ m/HIGH_SINGLE_MP/) {$splitline[2] = 12};
		if ($line =~ m/HIGH_SPAN_MP/) {$splitline[2] = 13};
		if ($line =~ m/STRECH_MP/) {$splitline[2] = 14};
		print OUTfrc "$splitline[0] $splitline[3] $splitline[4] $splitline[2]\n";
	}
	close INfrc;
	close OUTfrc;
}


########################### Print stats#################################

print "Total number of reads\n";
print "$readcount\n";
print "Good normal PE connections (%relative)\n";
print $readcount/2-$countcends-$countcsamewl-$countcdiff,"(100%)\n";
print "Good end connections:\n";
print "$countcends(",sprintf("%.3f",($countcends/(2*$numcontigs*$enddist))/(($readcount/2-$countcends-$countcsamewl-$countcdiff)/$totallength)*100),"%)\n";
print "Bad connection different contigs:\n";
print "$countcdiff(",sprintf("%.3f",($countcdiff/($totallength))/(($readcount-$countcends/2-$countcsamewl-$countcdiff)/$totallength)*100),"%)\n";
print "Bad connection same contig:\n";
print "$countcsamewl(",sprintf("%.3f",($countcsamewl/($totallength))/(($readcount-$countcends/2-$countcsamewl-$countcdiff)/$totallength)*100),"%)\n";

close OUTcends;
close OUTcdiff;
close OUTcsamewl;
close OUTckary;
close OUTccov;	
close OUTccovE;	
close OUTccovD;	
close OUTccovS;	
close OUTcgc;	
close OUTcrules;	

exit;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "samfile|i:s","fastafile|f:s", "enddist|e:s", "minlength|m:s", "avgreadlength|a:s", "headersplit|s:s", "circoscbin|b:s", "minpesplit|p:s", "frc|g:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );
    
	#if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{'infile'} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

    circosviz.pl

=head1 COPYRIGHT

   copyright (C) 2012 Mads Albertsen

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

	Indicates connections between contigs/scaffolds and produces circos 
	files for visualization

=head1 SYNOPSIS

circosviz.pl  -i -f [-e -m -a -s -b -p -g] : version 1.0

 [-help -h]           Displays this basic usage information
 [-samfile -i]        SAM formated mapping file
 [-fastafile -f]      Fasta file of sequences used in the mapping
 [-enddist -e]        Reads must hit within this distance of the end to be considered an end hit (default: 500)
 [-minlength -m]      The minimum lenght of a contig to infer closed circle (default: 3000)
 [-avgreadlength -a]  The average readlength of the reads used to estimate read position (default: 126)
 [-headersplit -s]    The symbol used to split the header to make the read 1 and 2 headers identical (default: _)
 [-circoscbin -b]     The windowlength used for the circos coverage file (default: 10000) 
 [-minpesplit -p]     Minimum distance between PE reads to be considered a split pe read (default: 1000)
 [-frc -g]            Add FRCbam Features using the Feature.gff output of FRCbam 	  
=cut
