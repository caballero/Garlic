#!/usr/bin/perl

=head1 NAME

benchmarkPClouds.pl

=head1 DESCRIPTION

check the TP, FP, TN, FN in RepeatMasker predicted regions.

=head1 USAGE

benchmarkRepeatMasker.pl FASTA REGION > OUT

=head1 EXAMPLES

benchmarkRepeatMasker.pl masked_sequence.fa sequence.out > out

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

$ARGV[1] or die "use benchmarkRepeatMasker.pl FASTA REGION [SHUFFLED?] > OUT\n";

my $fasta    = shift @ARGV;
my $region   = shift @ARGV;
my $shuf     = shift @ARGV;
my $mask_seq = '';
my $pred_seq = '';
my $name_seq = $fasta; $name_seq  =~ s/\.fa$//;
my $fn       =  0; # False negatives
my $fp       =  0; # False positives
my $tn       =  0; # True negatives
my $tp       =  0; # True positives
my $nl       =  0; # Errors
my $b        =  1;
my $sn       = 'NA';
my $sp       = 'NA';
my $acc      = 'NA';
my $fdr      = 'NA';


open F, "$fasta" or die "cannot open file $fasta\n";
while (<F>) {
    chomp;
    unless (m/>/) {
		s/[\s\r\n]//g;
        $mask_seq .= $_;
    }
}
close F;

$mask_seq =~ s/bad/NNN/ig;
$mask_seq =  uc($mask_seq) if (defined $shuf);
$pred_seq =  uc($mask_seq);
my $len   = length $mask_seq;
my $rep   = $mask_seq =~ tr/acgt/acgt/;
my $repf  = sprintf ("%.4f", $rep / $len);

open R, "$region" or die "cannot open file $region\n";
while (<R>) {
    chomp;
    s/^\s+//;
    next unless (m/\d+\s+\d+/);
    my @line = split (/\s+/, $_);
    my $ini  = $line[5];
    my $end  = $line[6];
    my $s = substr ($pred_seq, $ini - $b, $end - $ini + $b);
    substr ($pred_seq, $ini - $b, $end - $ini + $b) = lc $s;
}
close R;

for (my $i = 0; $i < length $mask_seq; $i++) {
     my $m = substr ($mask_seq, $i, 1);
     my $p = substr ($pred_seq, $i, 1);
     if    ($m =~ m/[ACGT]/ and $p =~ m/[ACGT]/) { $tn++; }
     elsif ($m =~ m/[ACGT]/ and $p =~ m/[acgt]/) { $fp++; }
     elsif ($m =~ m/[acgt]/ and $p =~ m/[ACGT]/) { $fn++; }
     elsif ($m =~ m/[acgt]/ and $p =~ m/[acgt]/) { $tp++; }
     else                                        { $nl++; }
}


$sp  = sprintf ("%.4f", $tp / ($tp + $fp)) if ( ($tp + $fp) > 0);
$sn  = sprintf ("%.4f", $tp / ($tp + $fn)) if ( ($tp + $fn) > 0);
$acc = sprintf ("%.4f", ($tp + $tn) / ($tp + $tn + $fp + $fn)) if ( ($tp + $tn + $fp + $fn) > 0);
$fdr = sprintf ("%.4f", $fp / ($fp + $tp)) if ( ($fp + $tp) > 0);

print join "\t", "SeqName", "Length", "\%Repeats", "TP", "TN", "FP", "FN", "Ns", "Spec", "Sens", "Acc", "FDR";
print "\n";
print join "\t", $name_seq, $len, $repf, $tp, $tn, $fp, $fn, $nl, $sp, $sn, $acc, $fdr;
print "\n";

