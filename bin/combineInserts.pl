#!/usr/bin/perl

use strict;
use warnings;

my %data;
my %suma;
my $gc;
opendir DIR, "." or die;
while (my $dir = readdir DIR) {
    next unless ($dir =~ m/^chr/);
    next unless (-e "$dir/$dir.inserts.W1000.data");
    warn "doing $dir\n";
    open FH, "$dir/$dir.inserts.W1000.data" or die;
    while (<FH>) {
        chomp;
        if (m/^#/) {
            $gc = $_;
        }
        else {
			my ($rep1, $rep2, $frq, $num) = split (/\t/, $_);
            $data{$gc}{$rep1}{$rep2} += $num;
            $suma{$gc}{$rep1} += $num;
        }
    }
    close FH;
}

warn "writing data\n";

foreach $gc (sort keys %data) {
    print "$gc\n";
	foreach my $rep1 (sort keys %{ $data{$gc} }) {
		my $sum = $suma{$gc}{$rep1};
		foreach my $rep2 (sort keys %{ $data{$gc}{$rep1} }) {
			my $num = $data{$gc}{$rep1}{$rep2};
			my $frq = sprintf("%.4f", $num /$sum);
			print join "\t", $rep1, $rep2, $frq, $num;
			print "\n";
		}
	}
}
warn "done.\n";
