#!/usr/bin/perl

=head1 NAME

sizeGenscan.pl

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
my ($id, $gid, $type, $dir, $tmp, $ini, $end, $len, $name, $gsize, $msize, $csize);
while (<>) {
    if (m/^Sequence (.+?) :/) {
        $name = $1;
    }
    next unless (m/^\s*\d+/);
    chomp;
    s/^\s+//;
    ($gid, $type, $dir, $ini, $end, $len) = split (/\s+/, $_);
    $gid =~ s/\.\d+$//;
    $id = "$name.$gid"; 
    
    if ($ini > $end) { #flip coordinates
        $tmp = $ini;
        $ini = $end;
        $end = $tmp;
    }
    
    $end++; # 0-based
    
    # genomic start and genomic end
    if (defined $seq{$id}{'ini'} ) {
        $seq{$id}{'ini'} = $ini if ($ini < $seq{$id}{'ini'});
    }
    else {
        $seq{$id}{'ini'} = $ini;
    }
    
    if (defined $seq{$id}{'end'} ) {
        $seq{$id}{'end'} = $end if ($end > $seq{$id}{'end'});
    }
    else {
        $seq{$id}{'end'} = $end;
    }
    
    next if ($type eq 'Prom');
    $seq{$id}{'mRNA'} += $len;
    
    next if ($type eq 'PlyA');
    $seq{$id}{'coding'} += $len;
}

foreach my $id (keys %seq) {
    $gsize = $seq{$id}{'end'} - $seq{$id}{'ini'};
    $msize = $seq{$id}{'mRNA'};
    $csize = $seq{$id}{'coding'};
    print "$id\t$gsize\t$msize\t$csize\n";
}
