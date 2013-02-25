#!/usr/bin/env perl
###############################################################################
#
#    calc.gc.pl
#
#	 Calculates gc content in fasta files.
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

my $inputfile;
my $outputfile;

$inputfile = &overrideDefault("inputfile.fasta",'inputfile');
$outputfile = &overrideDefault("outputfile.txt",'outputfile');
 
my $line;
my $seq2;
my $dummy = 0;
my @array;
my %counts;

######################################################################
# CODE HERE
######################################################################

	
open(IN, $inputfile) or die;
open(OUT, ">$outputfile") or die;
print OUT "contig\tgc\n";

while (my $line = <IN>)  {
	if ($line =~ m/>/) {
		chomp $line;
		if ($dummy == 1){
			push (@array, "$seq2");	
		}
		$seq2 = "$line\t";
		$dummy =1;
	}
	else {
		chomp $line;
		$seq2 = $seq2.$line;
		}
}
push (@array, "$seq2");	          #to catch the last sequence

foreach my $sequence (@array){
	$counts{G} = 0;
	$counts{C} = 0;
	$counts{A} = 0;
	$counts{T} = 0;
	my @tseq = split("\t", $sequence);
	my @seq = split("", $tseq[1]);
	foreach my $nucleotide (@seq) {
			$counts{$nucleotide}++;
		}
	my $gc = ($counts{G}+$counts{C})/($counts{G}+$counts{C}+$counts{A}+$counts{T})*100;
	$tseq[0] =~ s/>//;
	print OUT "$tseq[0]\t",sprintf("%.2f",$gc),"\n";
}

close IN;
close OUT;
exit;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inputfile|i:s", "outputfile|o:s");
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

    calc.gc.pl

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

	Calculates gc content in fastafiles.

=head1 SYNOPSIS

script.pl  -i -o [-h]

 [-help -h]           Displays this basic usage information
 [-inputfile -i]      Input fasta file. 
 [-outputfile -o]     Outputfile.
 
=cut