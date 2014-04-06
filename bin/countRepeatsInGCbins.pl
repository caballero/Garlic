#!/usr/bin/perl

my %data;
my @gc;
my %rep = ('SIMPLE' => 0);
my ($rep, $cls, $gc, $cnt);
while (<>) {
    next if (m/^\n/);
    if (m/^#GC=(\d+-\d+)/) {
        $gc = $1;
        push @gc, $gc;
    }
    elsif (m/^SIMPLE/) {
        $data{$gc}{'SIMPLE'}++;
    }
    elsif (m/(.+?):(.+?):/) {
        $rep = $1;
        $cls = $2;
        $data{$gc}{$cls}++;
        $rep{$cls}++;
    }
}

print join "\t", "REPEAT              ", @gc;
print "\n";
foreach $rep (sort keys %rep) {
    print "$rep";
    my $spaces = 20 - (length $rep);
    print " " x $spaces;
    foreach $gc (@gc) {
        $cnt = 0;
        $cnt = $data{$gc}{$rep} if (defined $data{$gc}{$rep});
        print "\t$cnt";
    }
    print "\n";
}
