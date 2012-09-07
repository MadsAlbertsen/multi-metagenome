#!/usr/bin/env perl
###############################################################################
#
#    img.matrix.pl
#
#	 Slices and dices raw IMG data. Makes a large matrix with absence or presence
#    of cats. Uses a genome file and a cat file to make the subset matrix.
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

my $incat;
my $catcolumn;
my $outputfile;
my $dir;
my $ingenomes;
my $columnsum;
my $category;
my $abundance;

$incat = &overrideDefault(".txt",'incat');
$ingenomes = &overrideDefault("genomes.txt",'ingenomes');
$outputfile = &overrideDefault("outputfile.txt",'outputfile');
$dir = &overrideDefault(".",'dir');
$columnsum = &overrideDefault(-1,'columnsum');
$category = &overrideDefault("pfam",'category');
$abundance = &overrideDefault(0,'abundance');
$catcolumn = &overrideDefault(9,'catcolumn');

my %cat;
my %genome;
my $outstr;
my %condensed;
my $genomecolumnsums;

######################################################################
# CODE HERE
######################################################################


open(INcat, $incat) or die("Cannot read file: $incat\n");
open(INgenomes, $ingenomes) or die("Cannot read file: $ingenomes\n");
open(OUT, ">$outputfile") or die("Cannot create file: $outputfile\n");

while ( my $line = <INgenomes> ) {                                                                 #Read in the genome file and store the information for subsequent parsing. Also defines which genomes are to be used.
	chomp $line;   
	my @splitline = split(/\t/,$line);
	$genomecolumnsums = scalar @splitline;
	$genome{$splitline[0]} = $line;
}

close INgenomes;

$outstr = $genome{"taxon_oid"};

while ( my $line = <INcat> ) {                                                                    #Read in the subset of PFAMs to use.
	chomp $line;   	
	$cat{$line} = 0;
}

foreach my $key (sort keys %cat){                                                                 #Print the PFAM id to the file
	$outstr = $outstr."\t".$key;
}

close INcat;

print OUT "$outstr\n";

opendir(DIR, $dir) or die "Cannot open dir: $dir!";

while ( my $filename = readdir(DIR)){                                                              #Read in all the individual raw IMG annotation PFAM files. One by one.
	if ($filename =~/.$category./){			
		my @id = split(/\./,$filename);                                                            #Extract the genome id from the filename
		if (exists $genome{$id[0]}){		
			open(INgenome, $filename) or die("Cannot read file: $filename\n");
			$outstr = $genome{$id[0]};
			while (my $line = <INgenome>)  {	                                                                   	
				chomp $line;
				my @splitline = split(/\t/,$line);
				if (exists $cat{$splitline[$catcolumn]}){                                                  
					if ($abundance == 0){
						$cat{$splitline[$catcolumn]} = 1;		                                           #Presence/absence	
					}
					else{
						if ($cat{$splitline[$catcolumn]} eq ""){
							$cat{$splitline[$catcolumn]} = 1;		
						}
						else{
							$cat{$splitline[$catcolumn]}++;                               #Abundance
						}
					}
				}
			}		
			foreach my $key (sort keys %cat){
				$outstr = $outstr."\t".$cat{$key};
				$cat{$key} = 0;
				}	
			print OUT "$outstr\n";
		}
	}
}

close OUT;

#If reduce flag is on then reduce by columnsum X.

if ($columnsum != -1){
	
	open(OUT, "$outputfile") or die("Cannot create file: $outputfile\n");
	open(OUT2, ">$outputfile".".columnsum.$columnsum.txt") or die("Cannot create file: $outputfile".".columnsum.$columnsum.txt\n");
	my $header = 0;
	while ( my $line = <OUT> ) {
		chomp $line;  
		$header++;
		my @splitline = split(/\t/,$line);
		my $entry = $splitline[$columnsum];		
		for (my $count = 1; $count <= $genomecolumnsums; $count++)  {                                 #Remove the descriptions from each columnsums so only the numbers remain.
			shift @splitline;		
		}
		my $count = -1;
		foreach my $element (@splitline)  {
			$count++;				
			if (exists $condensed{$entry}{$count}){
				$condensed{$entry}{$count} = $condensed{$entry}{$count}+$element;	
			}
			else{
				$condensed{$entry}{$count} = $element;
			}		
		}			
	}
	foreach my $key1 ( keys %condensed ) {
		my $outstr = $key1;
		foreach my $key2 ( keys %{$condensed{$key1}} ) {
			$outstr = $outstr."\t$condensed{$key1}{$key2}";
		}
		print OUT2 "$outstr\n";
	}
}

close OUT;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "incat|n:s", "outputfile|o:s","ingenomes|g:s","columnsum|c:s","category|t:s","dir|d:s","category|t:s","abundance|a:s","catcolumn|l:s");
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

img.matrix.pl  -n -g [-a -c -d -o -t -l -h]

 [-help -h]           Displays this basic usage information
 [-incat -n]          File containing e.g. PFAM names to extract from the annotation files. 
 [-category -t]       Category file identifier. E.g. pfam or cog (default: pfam).
 [-catcolumn -l]      Column in the annotation files containing the category id (default: 9).
 [-ingenomes -g]      List of genomes (taxon_oid) to extract data from.
 [-columnsum -c]      Make an condensed output using columnsum X in the genome file (default: none).
 [-dir -d]            Location of the individual annotation files (default: .).
 [-outputfile -o]     Outputfile (default: outputfile.txt). 
 [-abundance -a]      Output the number of each pfam (1) instead of presence/absence (0) (default: 0). 
 
=cut