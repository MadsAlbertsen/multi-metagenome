#!/usr/bin/env perl
###############################################################################
#
#    extract.long.hits.from.blast.pl - version 1.0
#
#	 Extracts the longest database hit given a blast report and the database.
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

my $inblast;
my $indb;
my $outputfile;
my $minlength;


$inblast = &overrideDefault("inblast.txt",'inblast');
$indb = &overrideDefault("indb.fa",'indb');
$outputfile = &overrideDefault("outputfile.txt",'outputfile');
$minlength = &overrideDefault("200",'minlength');
 
my %contig; 
 
######################################################################
# CODE HERE
######################################################################


open(INblast, $inblast) or die("Cannot read file: $inblast\n");

while ( my $line = <INblast> ) {
	chomp $line;   	
	my @splitline = split(/\t/, $line);	
	if (exists($contig{$splitline[1]})){
		my @splitline1 = split(/\t/, $contig{$splitline[1]});
		my $oldlength = abs($splitline1[0]-$splitline1[1]);
		my $newlength = abs($splitline[-3]-$splitline[-4]);
		if ($newlength > $oldlength){
			$contig{$splitline[1]} = "$splitline[-3]\t$splitline[-4]";
		}
	}
	else{
		$contig{$splitline[1]} = "$splitline[-3]\t$splitline[-4]";
	}
}

close INblast;

open(INdb, $indb) or die("Cannot read file: $indb\n");
open(OUT, ">$outputfile") or die("Cannot create file: $outputfile\n");

my $count = 0;
my $outstr = "ERROR";
my $prevheader = "ERROR";
my $header = "ERROR";
my $seq;


while ( my $line = <INdb> ) {
	chomp $line; 
	if ($line =~ m/>/) {
		$prevheader = $header;
		$header = $line;
		my @splitline = split(/>/, $prevheader);
		if($count > 0){
			if (exists($contig{$splitline[1]})){
				my @splitline1 = split(/\t/, $contig{$splitline[1]});
				if ($splitline1[0] > $splitline1[1]){					
					$outstr = substr($seq,$splitline1[1],$splitline1[0]-$splitline1[1]);
					
				}
				else{
					$outstr = substr($seq,$splitline1[0],$splitline1[1]-$splitline1[0]);
				}				
				my $len = length($outstr);
				my $clen = length($seq);
				if ($len > $minlength-1){
					print OUT "$prevheader.$splitline1[0].$splitline1[1].$len.$clen\n";
					print OUT "$outstr\n";
				}
			}
		}
		$seq = "";
		$count++;
	}
	else{
		$seq = $seq.$line;
	}
}
my @splitline = split(/>/, $header);
if (exists($contig{$splitline[1]})){
	my @splitline1 = split(/\t/, $contig{$splitline[1]});
	if ($splitline1[0] > $splitline1[1]){					
		$outstr = substr($seq,$splitline1[1],$splitline1[0]-$splitline1[1]);
	}
	else{
		$outstr = substr($seq,$splitline1[0],$splitline1[1]-$splitline1[0]);
	}
	my $len = length($outstr);
	if ($len > $minlength-1){	
		print OUT "$header.$splitline1[0].$splitline1[1].$len\n";
		print OUT "$outstr\n";
	}
}


close OUT;
close INdb;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inblast|b:s", "indb|d:s", "outputfile|o:s", "minlength|m:s");
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

    vprobes.generateprobes.pl

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

extract.long.hits.from.blast.pl  -i [-h] - version 1.0

 [-help -h]           Displays this basic usage information
 [-inblast -b]        Input blast file (outfmt 6).
 [-indb -d]           Input fastafile of the database used.
 [-outputfile -o]     Outputfile. 
 [-minlength -m]      Min. length to report sequence (default: 200).
 
=cut