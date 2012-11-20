#!/usr/bin/perl

=head1 NAME

runPclouds.pl

=head1 DESCRIPTION

Wrapper to run automatically PClouds program.

=head1 USAGE

perl runPclouds.pl -i FASTA -c CLOUDS -k KMER -o OUTPUT

=head1 EXAMPLES

perl runPclouds.pl -i test.fa -c hg19.w16_C10 -k 16 -o test.regions

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

=head1 CONTACT

jcaballero@systemsbiology.org

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with code.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Default parameters
my $help        = undef;         # Print help
my $verbose     = undef;         # Verbose mode
my $version     = undef;         # Version call flag
my $fasta       = undef;
my $clouds      = undef;
my $out         = undef;
my $kmer        = undef;

# Main variables
my $our_version = 0.1;        # Script version number
my $pclouds     = '/proj/hoodlab/juan/FakeSequence/PClouds/pcloud';

# Calling options
GetOptions(
    'h|help'            => \$help,
    'v|verbose'         => \$verbose,
    'version'           => \$version,
    'i=s'               => \$fasta,
    'c=s'               => \$clouds,
    'o=s'               => \$out,
    'k=i'               => \$kmer
) or pod2usage(-verbose => 2);

printVersion() if (defined $version);
pod2usage(-verbose => 2) if  (defined $help);
pod2usage(-verbose => 2) if !(defined $fasta);
pod2usage(-verbose => 2) if !(defined $clouds);
pod2usage(-verbose => 2) if !(defined $out);

open  C, ">Controlfile" or die "cannot write Controlfile\n";
print C <<__FILE__
Word size:
#OligoSize -> $kmer

Parameters to build P clouds:
#COPYTHRESHOLD -> 2
#ENDTHRESHOLD -> 10
#STEP1THRESHOLD -> 20
#STEP2THRESHOLD -> 200
#STEP3THRESHOLD -> 2000
#CALCHUNCKSIZE -> 10000000
#GENOMESIZE -> 3160000000

Parameters to define the repeat region:
#WindowSize -> 10
#PercentCutoff -> 80

Inputfile:
Input file to calculate counts:
#CountGenome -> hg19.pre

Input file to build p clouds and dissection
#OligoSets -> words$kmer.txt
#GenomeInput -> $fasta

Outputfile:
Clouds associated:
#MaincloudsInfo -> $clouds.mainclouds.info
#MaincloudsAssign -> $clouds.mainclouds.assign
#AcccloudsAssign -> $clouds.accclouds.assign

Genome related:
#CloudAnnotation -> $out.cloud
#RepeatRegion -> $out.region

Selection of Program parameter: (0 = not run, 1 = run)
#CALCOUNTS -> 0
#GETPCLOUDS -> 0
#DISSECTION -> 1
__FILE__
;
close C;

system ("$pclouds");

###################################
####   S U B R O U T I N E S   ####
###################################

# printVersion => return version number
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

