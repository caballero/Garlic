#!/usr/bin/perl

=head1 NAME

searchIntergenicBlocks.pl

=head1 DESCRIPTION

Search in a hard masked fasta file (N=N, X=genes, R=interspersed repeats, 
S=simple repeats) for regions with a specific length.

Output report the sequence name, coordinates (1-based) and GC content.

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -f --fasta       Input                      File*      STDIN
    -s --size        Fragment size              INT**      10kb 
    -o --output      Output                     File       STDOUT
    -h --help        Print this screen
    -v --verbose     Verbose mode
    
     * file can be compressed with gzip/bzip2
    ** kb, Mb and Gb supported 

=head1 EXAMPLES

    perl searchIntergenicBlocks.pl < FASTA > BLOCKS

    perl searchIntergenicBlocks.pl -f FASTA -o BLOCKS -s 1Mb

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
my $help       =  undef;
my $verbose    =  undef;
my $fasta      =  undef;
my $output     =  undef;
my $size       = '10kb';

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'f|fasta:s'       => \$fasta,
    's|size:s'        => \$size,
    'o|output:s'      => \$output
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $fasta) {
    my $fastah = $fasta;
    $fastah = "gunzip  -c $fasta | " if ($fasta =~ m/gz$/);    
    $fastah = "bunzip2 -c $fasta | " if ($fasta =~ m/bz2$/);
    open STDIN, "$fastah" or die "Cannot read file $fasta\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

my $fac = 1;
if ($size =~ m/k/i) {
    $fac = 1e3;
}
elsif ($size =~ m/m/i) {
    $fac = 1e6;
}
elsif ($size =~ m/g/i) {
    $fac = 1e9;
}
$size =~ s/\D//g;
$size *= $fac;
my $win = int($size / 10);

warn "loading sequences\n";
my %seq = ();
my ($id, $seq); 
while (<>) {
    chomp;
    if (m/^>/) {
        s/>//;
        $id = $_;
    }
    else {
        s/N/X/g;
        s/S/R/g;
        $seq{$id} .= $_;
    }
}

warn "searching regions\n";
foreach $id (keys %seq) {
    warn "  processing $id\n";
    my @frag = split (/(X+)/, $seq{$id});
    my $ini = 0;
    foreach my $frag (@frag) {
        my $len = length $frag;
        unless ($frag =~ m/X/) {
            if ($len >= $size) {
                my $n = $frag =~ tr/ACGTacgt/ACGTacgt/;
                next if (($n / $len) < 0.4);
                my $end = $ini + $size;
                print "$id\t$ini\t$end\\n";
            }
        }
        $ini += $len;
    }
}


