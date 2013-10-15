#!/usr/bin/env perl
###############################################################################
#
#    calc.coverage.in.bam.depth.pl
#    
#    Calculates coverage in depth file. The depth file can be generated using
#    bowtie2 and samtools.
#
#    bowtie2 -x assembly -U reads.fastq -S allignment.sam -p 10
#    bowtie2-build assembly.fa assembly
#    samtools view -bS allignment.sam > allignment.bam
#    samtools sort allignment.bam allignment.sorted.bam
#    samtools depth allignment.sorted.bam > depth.txt
#
#
#    Copyright (C) 2013 Mads Albertsen
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

my $outputfile;
my $indepth;

$indepth = &overrideDefault("depth.txt",'indepth');
$outputfile = &overrideDefault("output.txt",'outputfile');

my %coverage;
my @order;
my %length;
 
######################################################################
# CODE HERE
######################################################################


open(INdepth, "$indepth") or die("Cannot read file: $indepth\n");

while ( my $line = <INdepth> ) {
	chomp $line;
	my @splitline = split(/\t/,$line);
	if (exists($coverage{$splitline[0]})){
		$coverage{$splitline[0]} = $coverage{$splitline[0]} + $splitline[2];
		$length{$splitline[0]} = $splitline[1];
	}
	else{
		$coverage{$splitline[0]} = $splitline[2];
		$length{$splitline[0]} = $splitline[1];
		push (@order , $splitline[0]);
	}	
}

close INdepth;

open(OUT, ">$outputfile") or die("Cannot create file: $outputfile\n");

print OUT "Name,Average.coverage,Reference.length\n";

foreach my $scaffold (@order){
	my $cov = sprintf("%.3f", $coverage{$scaffold} / $length{$scaffold});
	print OUT "$scaffold,$cov,$length{$scaffold}\n";
}

close OUT;


######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "indepth|i:s", "outputfile|o:s");
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

    calc.coverage.in.bam.depth.pl

=head1 COPYRIGHT

   copyright (C) 2013 Mads Albertsen

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

script.pl  -i [-h]

 [-help -h]           Displays this basic usage information
 [-indepth -i]        Depth file of coverage generated using "samtools depth"
 [-outputfile -o]     Outputfile.
 
=cut
