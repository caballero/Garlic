#!/usr/bin/perl

use strict;
use warnings;

my %data;
my %suma;
my $gc;
my %total;
my $kmer = 'K4';

opendir DIR, "." or die;
while (my $dir = readdir DIR) {
    next unless ($dir =~ m/^chr/);
    next unless (-e "$dir/$dir.kmer.$kmer.W1000.data");
    warn "doing $dir\n";
    open FH, "$dir/$dir.kmer.$kmer.W1000.data" or die;
    while (<FH>) {
        chomp;
        if (m/^#/) {
            $gc = $_;
        }
        else {
			my ($word, $frq1, $frq2, $num) = split (/\s+/, $_);
			my $w = $word;
            my $b = chop $w;
            $data{$gc}{$w}{$b} += $num;
            $suma{$gc}{$w} += $num;
            $total{$gc} += $num;
        }
    }
    close FH;
}

warn "writing data\n";

foreach $gc (sort keys %data) {
    print "$gc\n";
	foreach my $w (sort keys %{ $data{$gc} }) {
		my $sum = $suma{$gc}{$w};
		foreach my $b (sort keys %{ $data{$gc}{$w} }) {
			my $num  = $data{$gc}{$w}{$b};
			my $frq1 = sprintf("%.8f", $num /$sum);
			my $frq2 = sprintf("%.8f", $num /$total{$gc});
			print join "\t", "$w$b", $frq1, $frq2, $num;
			print "\n";
		}
	}
}
warn "done.\n";
