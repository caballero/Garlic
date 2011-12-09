#!/usr/bin/perl

=head1 NAME

maskGenes.pl

=head1 DESCRIPTION


=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input                      File       STDIN
    -o --output      Output                     File       STDOUT
    -g --genes       Gene table                 File       genes.block
    -h --help        Print this screen
    -v --verbose     Verbose mode

=head1 EXAMPLES


=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2011

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
use Getopt::Long;
use Pod::Usage;

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $input      = undef;
my $output     = undef;
my $genes      = 'genes.block';

# Global variables
my %block = ();
my %seq   = ();
my ($seq, $seq_id);

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'g|genes:s'       => \$genes,
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $input) {
    open STDIN, "$input" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

warn "loading gene coordinates\n" if (defined $verbose);
open G, "$genes" or die "cannot open $genes\n";
while (<G>) {
    chomp;
    my ($seq_id, $ini, $end) = split (/\t/, $_);
    $block{$seq_id}{$ini}    = $end - $ini;
}
close G;

warn "loading sequences\n" if (defined $verbose);
while (<>) {
    chomp;
    if (m/>/) {
        s/>//;
        $seq_id = $_;
    }
    else {
        $seq{$seq_id} .= $_;
    }
}

while (($seq_id, $seq) = each %seq) {
    warn "masking regions in $seq_id\n" if (defined $verbose);
    if (defined $block{$seq_id}) {
        foreach my $pos (keys %{ $block{$seq_id} }) {
            my $len = $block{$seq_id}{$pos};
            substr($seq, $pos, $len) = 'N' x $len;
        }
    }
    
    $seq =~ tr/acgt/NNNN/; #repeats
    
    warn "writing sequence\n" if (defined $verbose);
    print ">$seq_id\n";
    while ($seq) {
        my $s = substr($seq, 0, 80);
        print "$s\n";
        substr($seq, 0, 80) = '';
    }
}
