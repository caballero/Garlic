#!/usr/bin/perl

=head1 NAME

sizeTwinscan.pl

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

=head1 CONTACT

jcaballero@systemsbiology.org

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with code.  If not, see <http://www.gnu.org/licenses/>.

=cut 

use strict;
use warnings;

my %seq = ();
my ($sid, $type, $prog, $ini, $end, $gsize, $csize);
while (<>) {
    next if (m/^#/);
    chomp;
    ($sid, $prog, $type, $ini, $end) = split (/\t/, $_);
    m/gene_id "(.+?)"/;
    $sid = $1;
    $sid =~ s/fasta.0+//;
    
    $end++;
    
    if (defined $seq{$sid}{'ini'}) {
       $seq{$sid}{'ini'} = $ini if ($ini < $seq{$sid}{'ini'});
    } 
    else {
       $seq{$sid}{'ini'} = $ini;
    }
    
    if (defined $seq{$sid}{'end'}) {
       $seq{$sid}{'end'} = $end if ($end > $seq{$sid}{'end'});
    } 
    else {
       $seq{$sid}{'end'} = $end;
    }
    
    $seq{$sid}{'coding'} += $end - $ini; 
}

foreach my $id (keys %seq) {
    $gsize = $seq{$id}{'end'} - $seq{$id}{'ini'};
    $csize = $seq{$id}{'coding'};
    print "$id\t$gsize\t$csize\n";
}
