#!/usr/bin/perl

=head1 NAME

getRegions.pl

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
use List::Util qw/shuffle/;

my $file = 'hg19.100kb.intergenic';
my $num  = 99;
my $dir  = './fasta';
my $len  = 1e5;
my %seq  = ();
my @reg  = ();
my ($seq, $reg, $ini, $end, $gc, $chr);

warn "loading sequences\n";
opendir D, "$dir" or die;
while (my $fa = readdir D) {
    next unless ($fa =~ m/\.fa$/);
    warn "    reading $fa\n";
    open F, "fasta/$fa" or die;
    while (<F>) {
        chomp;
        if (m/^>/) {
            s/>//;
            $chr = $_;
        }
        else {
            $seq{$chr} .= $_;
        }
    }
    close F;
}
closedir D;

warn "getting candidates regions\n";
open F, "$file" or die "cannot open $file\n";
@reg = <F>;
close F;
chomp(@reg);
@reg = shuffle(@reg);

warn "writing sequences\n";
for (my $i = 0; $i <= $num; $i++) {
    $reg = shift @reg;
    ($chr, $ini, $end, $gc) = split $reg;
    next unless (defined $seq{$chr});
    $seq = substr ($seq{$chr}, $ini - 1, $len);
    $end = $ini + $len;
    my $x = $i;
    $x = "0$i" if ($i < 10);
    warn "    intergenic_$x\n";
    open  F, "intergenic/>intergenic_$x.fa" or die "cannot open intergenic_$x.fa\n";
    print F ">intergenic_$x $chr:$ini-$end\n";
    while ($seq) {
        print F substr ($seq, 0, 70), "\n";
        substr ($seq, 0, 70) = '';
    }
    close F;
}

warn "done\n";
