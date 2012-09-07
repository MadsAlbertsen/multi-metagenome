#!/usr/bin/env perl
###############################################################################
#
#    cytoscape.extract.sub.graph.using.list.pl
#
#	 Given a list of nodes extracts all parts of the relating graph in a 
#    cytoscape connection file (nodes in column 0 and 2).
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

my $inconnections;
my $inlist;

$inconnections = &overrideDefault("incon.txt",'inconnections');
$inlist = &overrideDefault("inlist.txt",'inlist');

my %contigs;

######################################################################
# CODE HERE
######################################################################


open(INlist, $inlist) or die("Cannot read file: $inlist\n");
open(INcon, $inconnections) or die("Cannot read file: $inconnections\n");
open(OUT, ">$inlist.sub.txt") or die("Cannot create file: $inlist.sub.txt\n");
open(OUTsub, ">$inconnections.sub.txt") or die("Cannot create file: $inconnections.sub.txt\n");
open(OUTorg, ">$inlist.orginal.paint.cyto.txt") or die("Cannot create file: $inlist.orginal.paint.cyto.txt\n");

print OUTsub "node1\tinteraction\tnode2\tconnections\n";
print OUTorg "OrgScaffolds\n";

while ( my $line = <INlist> ) {
	chomp $line;   	
	$contigs{$line} = 1;
	print OUTorg "$line = 1\n";
}

close INlist;
close OUTorg;

my $newfound = 1;
my $count = 0;

while ($newfound == 1){
	$newfound = 0;
	$count++;
	print "Pass $count\n";
	while ( my $line = <INcon> ) {
		chomp $line;		
		my @splitline = split("\t",$line);
		if (exists($contigs{$splitline[0]}) and !exists($contigs{$splitline[2]}) ){
				$contigs{$splitline[2]} = 1;							
				$newfound = 1;
			}
		if (exists($contigs{$splitline[2]}) and !exists($contigs{$splitline[0]}) ){
				$contigs{$splitline[0]} = 1;							
				$newfound = 1;
			}
	}
	 seek(INcon,0,0);
}

foreach my $key (keys %contigs){
	print OUT "$key\n";
}

while ( my $line = <INcon> ) {
	chomp $line;		
	my @splitline = split("\t",$line);
	if (exists($contigs{$splitline[0]})){
		print OUTsub "$line\n";
	}
}

close INcon;
close OUT;
close OUTsub;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inlist|l:s", "inconnections|c:s");
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

script.pl  -i [-h]

 [-help -h]           Displays this basic usage information
 [-inlist -l]         List of nodes in subgraph to extract.
 [-inconnections -c]  Cytoscape connection file.
 
=cut