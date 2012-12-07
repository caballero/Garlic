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
    -b --base        Write base sequence here   File
    -r --repeat      write repeats info here    File
    -d --dir         Directory of models        Dir        ./data
    -w --win         Window size                INT        1000
    -k --kmer        Kmer size                  INT        8
    -b --block       Block size                 INT        10kb
    -m --model       Model to use*              STR        hg19
	-n --name        Temp files name            STR        fake
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
my $help       =   undef;
my $verbose    =   undef;
my $input      =   undef;
my $output     =   undef;
my $kmer       =       8;
my $win        =    1000;
my $block      =  '10kb';
my $model      =  'hg19';
my $dir        =  'data';
my $base       =   undef;
my $repeat     =   undef;
my $name       =  'fake';

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'k|kmer:i'        => \$kmer,
    'w|win:i'         => \$win,
    'b|block:s'       => \$block,
    'm|model:s'       => \$model,
    'a|base:s'        => \$base,
    'd|dir:s'         => \$dir,
    'r|repeat:s'      => \$repeat,
	'n|name:s'        => \$name
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $input) {
    my $infh = $input;
	$infh = "gunzip  -c $input | " if ($input =~ m/.gz$/);
	$infh = "bunzip2 -c $input | " if ($input =~ m/.bz2$/);
    open STDIN, "$infh" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

# Path to the createFakeSequence.pl script
my $creator = 'bin/createFakeSequenceBin10.pl';

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
    my $new_seq   = undef;
    my $bas_seq   = undef;
    my $rep_info  = undef;
    my $total_len = length $seq;
    for (my $i = 0; $i <= $total_len - $block; $i += $block) {
        my $block_size = $block;
        if ( ($i + $block) > $total_len) {
            $block_size = $total_len - (length $new_seq);
        }
        
        my $slice = substr ($seq, $i, $block_size);
        
        if ($slice =~ m/^N+$/) {
            $new_seq .= 'N' x $block;
            $bas_seq .= 'N' x $block;
            next;
        }
        
        my ($min, $max) = minmaxGC($slice);
        my $rep         = calcRepBases($slice) - 2;
        my $sim         = 2;
        $rep = 95 if ($rep > 95);
        $rep = 0  if ($rep <  1);
        warn "creating subseq $i GC=($min, $max), RF=$rep\n" if (defined $verbose);
        system ("perl $creator -m $model -k $kmer -w $win -r $rep -l $sim -g $min -c $max -s $block --write_base -n $name -d $dir");
        if (-e "$name.fasta") {
            my $new = '';
            open F, "$name.fasta" or die "cannot open $name.fasta\n";
            while (<F>) {
                chomp;
                next if (m/>/);
                $new .= $_;
            }
            close F;
            if ((length $new) > 1) {
                $new_seq .= $new;
                if (defined $base) {
                    open B, "$name.base.fasta" or die "cannot open $name.base.fasta\n";
                    while (<B>) {
                        next if (m/>/);
                        chomp;
                        $bas_seq .= $_;
                    }
                    close B;
                }
                if (defined $repeat) {
                    open R, "$name.inserts" or die "cannot open $name.inserts\n";
                    while (<R>) {
                        next if (m/POS/);
                        chomp;
                        my ($pos, $info) = split (/\t/, $_);
                        $pos+= $i;
                        $rep_info .= "$pos\t$info\n";
                    }
                    close R;
                }
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
    
    warn "masking N regions\n" if (defined $verbose);
    if ($seq =~ m/(N+)/g) {
        my $m = 0;
        foreach my $pos (@-) {
            my $nnum = $+[$m] - $-[$m];
            substr($new_seq, $pos, $nnum) = 'N' x $nnum;
            substr($bas_seq, $pos, $nnum) = 'N' x $nnum;
        }
    }
    
    # trimming if longer than the expected sequence
    $new_seq = substr($new_seq, 0, $total_len);
    $bas_seq = substr($bas_seq, 0, $total_len);
    
    warn "writing final sequence\n" if (defined $verbose);
    while ($new_seq) {
        print substr ($new_seq, 0, 70), "\n";
        substr ($new_seq, 0, 70) = '';
    }
    if (defined $base) {
        open  B, ">$base" or die "cannot open $base\n";
        print B ">$id (base) synthetic\n";
        while ($bas_seq) {
            print B substr ($bas_seq, 0, 70), "\n";
            substr ($bas_seq, 0, 70) = '';
        }
        close B;
    }
    if (defined $repeat) {
        open  R, ">$repeat" or die "cannot open $repeat\n";
        print R "POS\tINFO\n$rep_info";
        close R;
    }
}

# cleaning up
unlink "$name.fasta";
unlink "$name.base.fasta";
unlink "$name.inserts";
            
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
    $min_gc =  10 if ($min_gc == 1000);
    $max_gc = 100 if ($max_gc ==   -1);
    return ($min_gc, $max_gc);
}

sub calcRepBases {
    my $seq = shift @_;
    my $rep = int(rand 100);
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
    my $pgc = undef;
    
    if ($sum >= $win / 10) { # at least 10% of the sequence is useful
        $pgc  = int(100 * $ngc / $sum);
    }
    else { # no bases, use a random value
        my @gc = (10, 20, 30, 40, 40, 40, 50, 50, 50, 60, 60, 60, 70, 80, 90, 100);
        $pgc   = $gc[int(rand @gc)];
    }
    
	return $pgc;
}

#END
