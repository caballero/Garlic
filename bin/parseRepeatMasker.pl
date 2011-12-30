#!/usr/bin/perl

#parseRepeatMasker.pl < RM.out > LENGTHS
# simple script to parse the RM output, extracting Alu's lengths.
# Juan Caballero @ ISB 2011

while (<>) {
    s/^\s+//;
	next unless (m/^\d/);
	next unless (m/Alu/);
	my @line = split (/\s+/, $_);
	my $ini  = $line[5];
	my $end  = $line[6];
	my $len  = $end - $ini;
	next unless ($len > 0);
	print "$len\n";
}