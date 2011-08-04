#!/usr/bin/perl

=head1 NAME

extractUltraConservedRegions.pl

=head1 DESCRIPTION

Read a UCSC's multiz46.txt file, and filter the most conservated regions.
Output is a tab-delimited text table with columns:
1. UCSC id number
2. Chromosome
3. Direction (always +);
4. Position ini
5. Position end
6. Conservation score

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input                      File       STDIN
    -o --output      Output                     File       STDOUT
    -s --score       Score threshold            Float      0.65
    -l --length      Minimal length             Int        100
    -h --help        Print this screen
    -v --verbose     Verbose mode

=head1 EXAMPLES

  zcat multiz46ways.txt.gz | perl extractUltraConservedRegions.pl > ucr.txt
  perl extractUltraConservedRegions.pl -i multiz46ways.txt.gz -o ucr.txt
  perl extractUltraConservedRegions.pl -i multiz46ways.txt.gz -o ucr.txt -s 0.8

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2011

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

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $input      = undef;
my $output     = undef;
my $score      = 0.65;
my $length     = 100;

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    's|score:s'       => \$score,
    'l|length:i'      => \$length
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $input) {
    my $fileh = $input;
    $fileh = "gunzip  -c $input | " if ($input =~ m/gz$/);
    $fileh = "bunzip2 -c $input | " if ($input =~ m/bz2$/); 
    open STDIN, "$fileh" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

while (<>) {
    chomp;
    my ($id, $chr, $ini, $end, @rest) = split (/\t/, $_);
    my $sco = pop @rest;
    next unless ($sco >= $score);
    next unless ($end - $ini >= $length);
    print join "\t", $id, $chr, '+', $ini, $end, $sco;
    print "\n";
}

