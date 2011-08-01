#!/usr/bin/perl

=head1 NAME

correctKmerFreq.pl

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input                      File       STDIN
    -o --output      Output                     File       STDOUT
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

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $input) {
    open STDIN, "$input" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

my ($gc, $p, $q);
while (<>) {
    chomp;
    if (m/#GC=\d+-(\d+)/) {
        $gc = $1;
        $p  = $gc - 5;
        $q  = 100 - $p;
        print "$_\n";
    }
    else {
        my ($kmer, $frq, $cnt) = split (/\t/, $_);
        if ($cnt == 0 and $frq eq '0.25000000') {
            if ($kmer =~ m/[AT]$/) {
                $frq = $q / 2;
            } else {
                $frq = $p / 2;
            }
            print "$kmer\t$frq\t$cnt\n";
        }
        else {
            print "$_\n";
        }
    }
}
