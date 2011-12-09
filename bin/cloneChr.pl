#!/usr/bin/perl

=head1 NAME

cloneChr.pl

=head1 DESCRIPTION


=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input                      File       STDIN
    -o --output      Output                     File       STDOUT
    -b --bin         Window bin size            INT        10kb
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
my $bin        = 10000;

# Global Variables
my %seq;
my ($i, $n, $seq_id, $new_seq, $len, $seq, $slice, $new, $gc);

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'b|bin:s'         => \$bin
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $input) {
    open STDIN, "$input" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

warn "loading sequence\n" if (defined $verbose);
$/ = "\n>";
while (<>) {
    s/>//g;
    my @lines     = split (/\n/, $_);
    $seq_id       = shift @lines;
    $seq{$seq_id} = join ("", @lines);
}

while (($seq_id, $seq) = each %seq) {
    warn "cloning $seq_id\n" if (defined $verbose);
    $new_seq = '';
    while ($seq) {
        $slice = substr($seq, 0, $bin);
        $gc    = calcGC($slice);
        $len   = length $slice;
        system("bin/createFakeSequenceBin.pl -m hg19 -s $len -n new -g $gc -c $gc --no_repeat");
        $new   = '';
        open F, "new.fasta" or die "cannot open new.fasta\n";
        while (<F>) {
            next if (m/>/);
            chomp;
            $new .= $_;
        }
        close F;
        # masking
        if ($len != length $new) {
            my $new_len = length $new;
            die "lengths are different! want: $len, got: $new_len\n";
        }
        
        for ($i = 0; $i <= $len; $i++) {
            $n = substr($slice, $i, 1);
            substr($new, $i, 1) = 'N' unless ($n =~ m/[ATGC]/);
        }
        
        # piling and removing
        $new_seq .= $new;
        substr($seq, 0, $bin) = '';
    }
    unlink "new.fasta";
    
    warn "writing new sequence\n" if (defined $verbose);
    print ">$seq_id\n";
    while ($new_seq) {
        $seq = substr($new_seq, 0, 80);
        print "$seq\n";
        substr($new_seq, 0, 80) = '';
    }
}
        
sub calcGC {
    my $seq    = shift @_;
    my $num_gc = $seq =~ tr/GCgc/GCgc/;
    my $num_at = $seq =~ tr/ATat/ATat/;
    my $total  = $num_gc + $num_at;
    my $gc     = 37;
    if ($total > 1) {
        my $fgc = $num_gc / $total;
        if    ($fgc <= 0.37) { $gc =  37; }
        elsif ($fgc <= 0.39) { $gc =  39; }
        elsif ($fgc <= 0.42) { $gc =  42; }
        elsif ($fgc <= 0.45) { $gc =  45; }
        else                 { $gc = 100; }
    }
    return $gc;
}
