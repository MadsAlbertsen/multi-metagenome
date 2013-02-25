#!/usr/bin/env perl
###############################################################################
#
#    extract.using.header.list.pl - version 1.0
#
#	 Extracts a subset of seauences given a list of the headers to extract.
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

my $inlist;
my $insequences;
my $outputfile;

$inlist = &overrideDefault("linlist.txt",'inlist');
$insequences = &overrideDefault("sequences.fa",'insequences');
$outputfile = &overrideDefault("outputfile.txt",'outputfile');
 
my $countp = 0;
my $countout = 0;
my $count = 0;
my $seq;
my $header = ">start[test]";
my $prevheader = ">start[test]";
my %extract;

 
######################################################################
# CODE HERE
######################################################################


open(INlist, $inlist) or die("Cannot read file: $inlist\n");

while ( my $line = <INlist> ) {
	chomp $line; 
		if (!exists($extract{$line})){
			$extract{$line} = 1;
			$count++;
		}		
}
close INlist;

print "$count sequences to extract.\n";


open(INseq, $insequences) or die("Cannot read file: $insequences\n");
open(OUT, ">$outputfile") or die("Cannot create file: $outputfile\n");

while ( my $line = <INseq> ) {
	chomp $line; 
	if ($line =~ m/>/) {
		$countp++;
		$prevheader = $header;
		$header = $line;
		my @splitline = split(/>/,$prevheader);		
		if(exists($extract{$splitline[1]}) or exists($extract{$prevheader})){
			print OUT "$prevheader\n";
			print OUT "$seq\n";
			$countout++;
		}
		$seq = "";
	}
	else{
		$seq = $seq.$line;
	}
}
my @splitline = split(/>/,$header);
if(exists($extract{$splitline[1]}) or exists($extract{$header})){
	print OUT "$header\n";
	print OUT "$seq\n";
	$countout++;
}

print "$countp sequences in total.\n";
print "$countout sequences extracted.\n";

close OUT;
close INseq;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "insequences|s:s", "inlist|l:s", "outputfile|o:s");
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

    extract.using.header.list.pl

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



=head1 SYNOPSIS

extract.using.header.list.pl  -i [-h] - version 1.0

 [-help -h]           Displays this basic usage information
 [-inlist -l]         List of headers to use for extraction.
 [-insequences -s]    Sequence file where a subset of sequences are to be extracted from.
 [-outputfile -o]     Outputfile. 
 
=cut