#!/usr/bin/perl

# parsePclouds.pl < REGION > LENGTHS
# simple script to parse Pclouds output (*region file)
# Juan Caballero @ ISB 2011

use strict;
use warnings;

while (<>) {
	chomp;
	my ($ini, $end) = split (/t/, $_);
	next unless (defined $end);
	my $len  = $end - $ini;
	next unless ($len > 0);
	print "$len\n";}