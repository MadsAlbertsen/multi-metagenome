#!/usr/bin/env perl
###############################################################################
#
#    extract.pe.reads.using.single.pl
#
#	 Given a list of single reads it extracts the PE reads from 2 PE fastq files.   
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


my $insingle;
my $inread;
my $splitheader;
my $splitfasta;
my $print = 0;
my $linecount = 0;
my $printcount = 0;

$inread = &overrideDefault("paired.fa",'inread');
$insingle = &overrideDefault("single.fa",'insingle');
$splitheader = &overrideDefault("_",'splitheader');
$splitfasta = &overrideDefault("_",'splitfasta');
 
my $seq = ""; 
my $readfound = 0;
my $toextract = 0;
my $extracted = 0;
my %reads;

######################################################################
# CODE HERE
######################################################################

open(INsingle, $insingle) or die("Cannot read file: $insingle\n");                                    #First read in all headers in the read 1 file that need to be matched in the read 2 file.

while ( my $line = <INsingle> ) {
	chomp $line;   	
	if ($line =~ m/>/) {
		my @splitline = split(/$splitfasta/,$line);		
		$reads{$splitline[0]} = 1;
		$toextract++;
	}
}
print "Found $toextract single reads.\n";
close INsingle;

open(OUT, ">paired.sub.fa") or die("Cannot create file: paired.sub.fa\n");
open(INread, $inread) or die("Cannot read file: $inread\n");

while (my $line = <INread>)  {	                                                                   #Look for matching read1 headers in the read2 file.
	chomp $line;
	if ($line =~ m/>/){
		$linecount++;
		$printcount++;
		my @splitline = split(/$splitheader/,$line);
		if (exists($reads{$splitline[0]})){
			$print = 1;
			$extracted++;
		}
		else{
		$print = 0;
		}
	}			
	if ($print == 1){
		print OUT "$line\n";
	}
	if ($printcount == 1000000){
		print "$linecount reads scanned - $extracted extracted\n";
		$printcount = 0;
	}
}
print "Extracted $extracted of $linecount reads.\n";
close INread;
close OUT;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inread|p:s","splitheader|x:s","insingle|s:s","splitfasta|y:s");
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

    extract.read2.using.read1.pl

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

Used in digital normalization. First read1 library is digital normalized
by khmer scripts and then read2 is extracted using the remaining read1 reads
using this scripts. This ensures proper use of PE reads. 


=head1 SYNOPSIS

extract.read2.using.read1.pl  -f -r -s [-h -x]

 [-help -h]           Displays this basic usage information
 [-inread -p]        pairedreads.fa.
 [-insingle -s]       Singlereads.fa. 
 [-splitheader -x]    Code used to split the header of the fastq files (default: "_")
 [-splitfasta -y]     Code used to split the header of the fasta file (default: " ")
 
=cut