#!/usr/bin/perl

=head1 NAME

kmerShuffle.pl

=head1 DESCRIPTION

Permute a fasta file with sequences based on kmer frequencies.

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input*                     File       STDIN
    -o --output      Output                     File       STDOUT
    -k --kmer        Kmer size                  Int        2
    -w --win         Window size**              Int        10kb
    -c --col         Column size                Int        80
    -r --repeat      Don't uppercase repeats
    -h --help        Print this screen
    -v --verbose     Verbose mode

     * Accepts compressed files with gzip/bzip2
    ** Can use kb, Mb, Gb symbols
    
=head1 EXAMPLES

    kmerShuffle.pl < FASTA > SHUFFLED

    kmerShuffle.pl -i FASTA.gz -o SHUFFLED -k 4 -w 1Mb
    
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
use List::Util 'shuffle';

# Parameters initialization
my $help       =  undef;
my $verbose    =  undef;
my $input      =  undef;
my $output     =  undef;
my $win        = '10kb';
my $kmer       =      2;
my $col        =     80;
my $repeat     =  undef;

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'w|window:s'      => \$win,
    'k|kmer:i'        => \$kmer,
    'c|col:i'         => \$col,
    'r|repeat'        => \$repeat
);

pod2usage(-verbose => 2) if (defined $help);

# File opening (if required)
if (defined $input) {
    my $inputh = $input;
    $inputh = "gunzip  -c $input | " if ($input =~ m/gz$/i);
    $inputh = "bunzip2 -c $input | " if ($input =~ m/bz2$/i); 
    open STDIN, "$inputh" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

# Window expansion
my $fac = 1;
if    ($win =~ m/k/i) { $fac = 1e3; }
elsif ($win =~ m/m/i) { $fac = 1e6; }
elsif ($win =~ m/g/i) { $fac = 1e9; }

$win =~ s/\D//g;
$win *= $fac;

# Kmer - 1 = Jmer
my $jmer = $kmer - 1;

# Reading fasta file
$/ = "\n>"; # slurp mode
while (<>) {
    s/>//g;
    my @seq  = split (/\n/, $_);
    my $name = shift @seq;
    warn "shuffling $name\n" if (defined $verbose);
    my $seq = join "", @seq;
    $seq =~ s/[^ACGTacgtNn]/N/g;
    $seq = uc $seq unless (defined $repeat);
    my $new = '';
    while ($seq) {
        my $s = substr ($seq, 0, $win);
        $new .= kshuffle($s);
        substr ($seq, 0, $win) = '';
    }
    print ">$name\n";
    while ($new) {
        print substr ($new, 0, $col), "\n";
        substr ($new, 0, $col) = '';
    }
}
warn "done\n" if (defined $verbose);

sub kshuffle {
    my $seq  = shift @_;
    my $len  = length $seq;
    my @kmer = ();
    my ($res, $num, $cnt, $nt, $ext, $seed);
    
    # kmer counting
    for (my $i = 0; $i <= $len - $kmer - 1; $i++) {
        push @kmer, substr ($seq, $i, $kmer);
    }

    # permute the kmers
    @kmer = shuffle(@kmer);
    
    # first element (random)
    $res  = shift @kmer;
    $seed = substr ($res, -$jmer);
    $num  = $#kmer;
    $cnt  = 0;
    while (1) {
        $ext = shift @kmer;
        # extend if we match the seed
        if ($ext =~ m/^$seed/) {
            $nt   = chop $ext;
            $res .= $nt;
            $seed = substr ($res, -$jmer);
            $cnt  = 0;
            $num--;
        }
        else {
            $cnt++;
            # continue if we reach all the possibles ways
            if ($cnt > $num) {
                $nt   = chop $ext;
                $res .= $nt;
                $seed = substr ($res, -$jmer);
                $cnt  = 0;
                $num--;
            }
            else {
                push @kmer, $ext;
            }
        }
        # stop when we finish all kmers
        last unless (defined $kmer[0]);
    }
    return $res;
}
