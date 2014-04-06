#!/usr/bin/perl 

use strict;
use warnings;

my %rep;
my %fam;
my $nseq = -1;
while (<>) {
	if (m/^#/) {
		$nseq++;
		next;
    }
	chomp;
	my ($ini, $end, $rep) = split (/\t/, $_);
	my @rep = split (/,/, $rep);
	foreach my $r (@rep) {
		if ($r =~ m/SIMPLE/) {
			$rep{$nseq}{'SIMPLE'}++;
			$fam{'SIMPLE'} = 1;
			
		}
		else {
			my ($type, $fam) = split (/:/, $r);
			$fam =~ s/\/.*//;
			$rep{$nseq}{$fam}++;
			$fam{$fam} = 1;
		}
	}
}

my @fam = sort keys %fam;
print join "\t", 'Seq', @fam;
print "\n";

for (my $i = 0; $i <= $nseq; $i++) {
	my $n = $i;
	$n = "0$n" if ($i < 10);
	print "fake$n";
	foreach my $fam (@fam) { 
		my $cnt = 0;
		$cnt = $rep{$i}{$fam} if (defined $rep{$i}{$fam});
		print "\t$cnt";
	}
	print "\n";
}
