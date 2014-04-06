#!/usr/bin/perl

use strict;
use warnings;

my %data;
my $gc;
opendir DIR, "." or die;
while (my $dir = readdir DIR) {
    next unless ($dir =~ m/^chr/);
    next unless (-e "$dir/$dir.repeats.W1000.data");
    warn "doing $dir\n";
    open FH, "$dir/$dir.repeats.W1000.data" or die;
    while (<FH>) {
        chomp;
        if (m/^#/) {
            $gc = $_;
        }
        elsif (m/:/) {
            push @{ $data{$gc} }, $_;
        }
    }
    close FH;
}

warn "writing data\n";

foreach $gc (sort keys %data) {
    print "$gc\n";
    print join "\n", @{ $data{$gc} };
    print "\n";
}
warn "done.\n";
