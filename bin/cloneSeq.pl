#!/usr/bin/perl

=head1 NAME

closeSeq.pl

=head1 DESCRIPTION

Read a fasta file and create a new sequence based in "de novo" fragments with
createFakeSequence.pl

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input Fasta                File       STDIN
    -o --output      Output Fasta               File       STDOUT
    -w --win         Window size                INT        1000
    -k --kmer        Kmer size                  INT        4
    -b --block       Block size                 INT        10kb
    -m --model       Model to use*              STR        hg19
    -h --help        Print this screen
    -v --verbose     Verbose mode

    * this variable uses the names defined by the UCSC genome browser: hg19
      mm9, fr2, ...
      
=head1 EXAMPLES


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
my $kmer       =     4;
my $win        =  1000;
my $block      = '10k';
my $model      = 'hg19';

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'k|kmer:i'        => \$kmer,
    'w|win:i'         => \$win,
    'b|block:s'       => \$block,
    'm|model:s'       => \$model
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $input) {
    open STDIN, "$input" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

# Path to the createFakeSequence.pl script
my $creator = 'bin/createFakeSequenceBin.pl';

# main varibles
my %seq = ();
my ($id, $seq, $gc_min, $gc_max);

# size adjustment
my $fac = 1;
if    ($block =~ m/k/i) { $fac = 1e3; }
elsif ($block =~ m/m/i) { $fac = 1e6; }
elsif ($block =~ m/g/i) { $fac = 1e9; } 
$block =~ s/\D//g;
$block *= $fac;

# processing sequences
warn "loading sequences\n" if (defined $verbose);
while (<>) {
    chomp;
    if (/>/) {
        s/>//;
        $id = $_;
    }
    else {
        $seq{$id} .= $_;
    }
}

warn "creating homolog sequences\n" if (defined $verbose);
while (($id, $seq) = each %seq) {
    print ">$id synthetic\n";
    my $new_seq = undef;
    for (my $i = 0; $i <= (length $seq) - $block; $i += $block) {
        my $slice = substr ($seq, $i, $block);
        my ($min, $max) = minmaxGC($slice);
        my $rep = calcRepBases($slice);
        warn "creating subseq $i GC=($min, $max), RF=$rep\n" if (defined $verbose);
        system ("perl $creator -m $model -k $kmer -w $win -r $rep -g $min -c $max -s $block -n fake");
        if (-e 'fake.fasta') {
            my $new = '';
            open F, "fake.fasta" or die "cannot open fake.fasta\n";
            while (<F>) {
                chomp;
                next if (m/>/);
                $new .= $_;
            }
            close F;
            #unlink 'fake.fasta';
            #unlink 'fake.inserts';
            if ((length $new) > 1) {
                $new_seq .= $new;
            }
            else {
                $i--;
                warn "failed (no sequence found), redoing\n" if (defined $verbose);
            }
        }
        else {
            $i--;
            warn "failed (no sequence found), redoing\n" if (defined $verbose);
        }
    }
    warn "writing final sequence\n" if (defined $verbose);
    while ($new_seq) {
        print substr ($new_seq, 0, 70), "\n";
        substr ($new_seq, 0, 70) = '';
    }
}

#################################################
##      S  U  B  R  O  U  T  I  N  E  S        ##
#################################################
sub minmaxGC {
    my $seq     = shift @_;
    my $min_gc  = 1000;
    my $max_gc  = -1;
    for (my $i = 0; $i <= (length $seq) - $win; $i += $win) {
        my $s = substr ($seq, $i, $win);
        my $gc = calcGC($s);
        $min_gc = $gc if ($gc < $min_gc);
        $max_gc = $gc if ($gc > $max_gc);
    }
    return ($min_gc, $max_gc);
}

sub calcRepBases {
    my $seq = shift @_;
    my $rep = rand;
    my $rrb = $seq =~ tr/acgt//;
    my $eeb = $seq =~ tr/ACGT//;
    my $sum = $rrb + $eeb;
    $rep = int(100 * $rrb / $sum) if ($sum > 1);
    return $rep;
}

sub calcGC {
    my $seq = shift @_;
    my $ngc = $seq =~ tr/GCgc//;
    my $nat = $seq =~ tr/ATat//;
    my $sum = $nat + $ngc;
    my $gc  = undef;
    my $pgc = undef;
    
    if ($sum >= $win / 10) { # at least 10% of the sequence is useful
        $gc = $ngc / $sum;
    }
    else { # no bases, use a random value
        my @gc = (0.3699, 0.3899, 0.4199, 0.4499, 1);
        $gc = $gc[int(rand @gc)];
    }
    
    if    ($gc <= 0.37) { $pgc =  37; }
    elsif ($gc <= 0.39) { $pgc =  39; }   
    elsif ($gc <= 0.42) { $pgc =  42; }
    elsif ($gc <= 0.45) { $pgc =  45; }
    else                { $pgc = 100; }
    
	return $pgc;
}

