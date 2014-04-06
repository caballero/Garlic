#!/usr/bin/perl

use strict;
use warnings;

my %data;
my %suma;
opendir DIR, "." or die;
while (my $dir = readdir DIR) {
    next unless ($dir =~ m/^chr/);
    next unless (-e "$dir/$dir.GCt.W1000.data");
    warn "doing $dir\n";
    open FH, "$dir/$dir.GCt.W1000.data" or die;
    while (<FH>) {
        chomp;
		my ($gc1, $gc2, $frq, $num) = split (/\t/, $_);
        $data{$gc1}{$gc2} += $num;
        $suma{$gc1} += $num;
    }
    close FH;
}

warn "writing data\n";

foreach my $gc1 (sort keys %data) {
	my $sum = $suma{$gc1};
	foreach my $gc2 (sort keys %{ $data{$gc1} }) {
		my $num = $data{$gc1}{$gc2};
		my $frq = sprintf("%.8f", $num /$sum);
		print join "\t", $gc1, $gc2, $frq, $num;
		print "\n";
	}
}
warn "done.\n";
