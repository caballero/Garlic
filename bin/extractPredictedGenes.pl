#!/usr/bin/perl

=head1 NAME

extractPredictedGenes.pl

=head1 DESCRIPTION

simple script to parse and extract the sequences from gene prediction programs
outputs (Augustus, Genscan, Twinscan)

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input                      File*      STDIN
    -o --output      Output                     File       STDOUT
    -p --program     Input type                 Prog**      
    -l --length      Print seq len only
    -b --basename    Use this basename          STR        sequence
    -n --num         Start counting here        INT        1
    -h --help        Print this screen
    -v --verbose     Verbose mode
    
     * File can be compresed with gz|bzip2
    ** Programs supported are (A)ugustus, (G)enscan and (T)winscan.

=head1 EXAMPLES

    extractPredictedGenes.pl -p A < augustus.out > predict.fa
    extractPredictedGenes.pl -i genscan.out -p G -l -o predict.sizes
    extractPredictedGenes.pl -i twinscan.out.gz -p T -o predict.fa

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
my $intype     = undef;
my $only_size  = undef;
my $base       = 'sequence';
my $num        = 0;

# Main variables
my ($program, $seq, $len, $rec);

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'p|program:s'     => \$intype,
    'l|length'        => \$only_size,
    'b|basename:s'    => \$base,
    'n|num:i'         => \$num
);
pod2usage(-verbose => 2) if (defined $help);

# Check the input format
$program = 'augustus' if ($intype =~ m/^a/i);
$program = 'genscan'  if ($intype =~ m/^g/i);
$program = 'twinscan' if ($intype =~ m/^t/i);

# Open files if required
if (defined $input) {
    open STDIN, "$input" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

# Run parsers
if    ($program eq 'augustus') { parseAugustus(); }
elsif ($program eq  'genscan') {  parseGenscan(); }
elsif ($program eq 'twinscan') { parseTwinscan(); }
else  {    die "Format not recognized $intype\n"; }

sub printFasta {
    my ($name, $seq) = @_;
    print ">$name\n";
    while ($seq) {
        my $s = substr($seq, 0, 80);
        print "$s\n";
        substr($seq, 0, 80) = '';
    }
}

sub parseAugustus {
    $seq = '';
    $len =  0;
    $rec =  0;
    while (<>) {
        if (m/# protein sequence =/) {
            $num++;
            $rec = 1;
            chomp;
            s/# protein sequence = \[//;
            s/\]//;
            $seq = $_;
        }
        elsif ($rec == 1) {
            if (m/### end gene/) {
                $rec = 0;
                $len = length $seq;
                if (defined $only_size) {
                    print "$len\n";
                }
                else {
                    printFasta("$base$num", $seq);
                }
                $seq = '';
            }
            else {
                s/# //;
                s/\]//;
                $seq .= $_;
            }
        }
        else {
            # do nothing
        }
    }
}

sub parseGenscan {
    $seq = '';
    $len =  0;
    $rec =  0;
    while (<>) {
        if (m/^>/) {
            $rec = 1;
            $num++;
        }
        elsif ($rec == 1) {
            if (m/^\n/) {
                $rec = 0;
                $len = length $seq;
                if (defined $only_size) {
                    print "$len\n";
                }
                else {
                    printFasta("$base$num", $seq);
                }
                $seq = '';
            }
            else {
                chomp;
                $seq .= $_;
            }
        }
        else {
            #do nothing
        }
    }
    # if the file have a last sequence
    if (defined $seq) {
        if (defined $only_size) {
            print "$len\n";
        }
        else {
            printFasta("$base$num", $seq);
        }
    }
}

sub parseTwinscan {
    warn "twinscan doesn't report the sequences, I can output the prediction sizes only\n";
    $len = 0;
    $rec = 'first_prediction';
    while (<>) {
        next if (m/^#/);
        my @line = split (/\t/, $_);
        next unless ($line[2] =~ m/CDS|start_codon|stop_codon/ );
        
        $line[-1] =~ m/gene_id "(.+?)"/;
        my $gen = $1;
        if ($gen eq $rec) {
            $len += $line[4] - $line[3];
        }
        else {
            print "$len\n" unless ($rec eq 'first_prediction');
            $len = 0;
            $rec = $gen;
        }
    }
    # last prediction
    print "$len\n" unless ($rec eq 'first_prediction');
}
    
