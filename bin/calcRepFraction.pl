#!/usr/bin/perl 

use strict;
use warnings;

$/ = "\n>";
while (<>) {
	s/>//g;
	my ($sid, @seq) = split (/\n/, $_);
	my $seq = join "", @seq;
	my $num_rep = $seq =~ tr/acgt/acgt/;
	print $num_rep / length $seq;
	print "\n";
}
