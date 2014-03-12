#!/usr/bin/perl 

use strict;
use warnings;

my %rep;
my %cnt;
my ($seq, $rep, $rid, $cnt);

while (<>) {
	s/^\s+//;
	chomp;
	next unless (m/^\d+/);
	my @line = split (/\s+/, $_);
	$seq  =  $line[4];
	$rep  =  $line[10]; 
	$rep  =~ s/\/.*//;
	$rep  =~ s/\?//;
	$rep  = "Simple_repeat" if ($rep eq "Low_complexity");
	$rid  =  $line[-1];
	$rep{$rep} = 1;
	$cnt{$seq}{$rep}{$rid}++;
}

my @rep = sort keys %rep;
print join "\t", "SEQ", @rep;
print "\n";

foreach $seq (sort keys %cnt) {
	print $seq;
	foreach $rep (@rep) {
		$cnt = 0;
		foreach $rid (keys %{ $cnt{$seq}{$rep} }) {
			$cnt++;
		}
		print "\t$cnt";
	}
	print "\n";
}
