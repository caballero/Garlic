#!/usr/bin/perl

=head1 NAME

benchmarkPClouds.pl

=head1 DESCRIPTION

check the TP, FP, TN, FN in PClouds predicted regions.

=head1 USAGE

benchmarkPClouds.pl FASTA REGION > OUT

=head1 EXAMPLES

benchmarkPClouds.pl masked_sequence.fa sequence.region > out

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

$ARGV[1] or die "use benchmarkPClouds.pl FASTA REGION > OUT\n";

my $fasta    = shift @ARGV;
my $region   = shift @ARGV;
my $mask_seq = '';
my $pred_seq = '';
my $name_seq = '';
my $fn       =  0; # False negatives
my $fp       =  0; # False positives
my $tn       =  0; # True negatives
my $tp       =  0; # True positives
my $nl       =  0; # Errors
my $b        =  0;

open F, "$fasta" or die "cannot open file $fasta\n";
while (<F>) {
    chomp;
    if (m/>(\S+)/) {
        $name_seq = $1; 
    }
    else {
        $mask_seq .= $_;
    }
}
close F;

$pred_seq = uc($mask_seq);

open R, "$region" or die "cannot open file $region\n";
while (<R>) {
    chomp;
    next unless (m/\d+\s+\d+/);
    my ($ini, $end) = split (/\s+/, $_);
    my $s = substr ($pred_seq, $ini - $b, $end - $ini);
    substr ($pred_seq, $ini - $b, $end - $ini) = lc $s;
}
close R;

for (my $i = 0; $i <= length $mask_seq; $i++) {
     my $m = substr ($mask_seq, $i, 1);
     my $p = substr ($pred_seq, $i, 1);
     if    ($m =~ m/[ACGT]/ and $p =~ m/[ACGT]/) { $tn++; }
     elsif ($m =~ m/[ACGT]/ and $p =~ m/[acgt]/) { $fp++; }
     elsif ($m =~ m/[acgt]/ and $p =~ m/[ACGT]/) { $fn++; }
     elsif ($m =~ m/[acgt]/ and $p =~ m/[acgt]/) { $tp++; }
     else                                        { $nl++; }
}

my $sp = sprintf ("%.4f", $tp / ($tp + $fp));
my $sn = sprintf ("%.4f", $tp / ($tp + $fn));

print join "\t", $name_seq, $tp, $tn, $fp, $fn, $nl, $sp, $sn;
print "\n";

