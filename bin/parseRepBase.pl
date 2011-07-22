#!/usr/bin/perl

=head1 NAME

parseRepBase.pl

=head1 DESCRIPTION

Search for RepBase fasta files (*.ref) and remove duplicates, creating a fasta 
file with unique sequences.

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -d --dir         Working directory          Dir        LOCAL
    -o --output      Output                     File       STDOUT
    -e --exclude     Files to exclude (RegEx)   Str        simple
    -h --help        Print this screen
    -v --verbose     Verbose mode

=head1 EXAMPLES

   perl parseRepBase.pl > RepBase.fa
   
   perl parseRepBase.pl -d /opt/RepBase -o RepBase.fa

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
use File::Find;

# Parameters initialization
my $help       =    undef;
my $verbose    =    undef;
my $dir        =      '.'; # local dir is default
my $output     =    undef;
my $exclude    = 'simple'; # no simple repeats

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'd|dir:s'         => \$dir,
    'o|output:s'      => \$output,
    'e|exclude:s'     => \$exclude
);
pod2usage(-verbose => 2) if (defined $help);

if ($dir eq '.') {
    warn "searching files in local directory\n" if (defined $verbose);
}
else {
    warn "searching files in directory $dir\n" if (defined $verbose);
}

my @files = split (/,/, searchFiles('.ref', $dir));
my %seq   = ();
$exclude  = s/,/|/g;

foreach my $file (@files) {
    next if ($file =~ m/$exclude/);
    warn "loading sequence in $file\n" if (defined $verbose);
    open F, "$file" or die "cannot read $file\n";
    my $id = undef;
    $/ = "\n>"; # slurp mode per sequence
    while (<F>) {
        chomp;
        s/>//;
        my @seq = split (/\n/, $_);
        my $name = shift @seq;
        $name =~ m/(.+?)\t(.+?)\t/;
        my $type = $1;
        my $fam  = $2;
        $id = "$type:$fam";
        next if (defined $seq{$id});
        $seq{$id} = join "\n", ">$id", @seq;
    }
    close F;
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

my $cnt = 0;
while ( my ($rep, $seq) = each %seq ) {
    print $seq;
    $cnt++;
}
warn "$cnt sequences found\n" if (defined $verbose);

sub searchFiles {
    my ($pat, $dir) = @_;
    my $files = undef;
    find ( sub { $files .= $File::Find::name . ',' if (m/$pat/) }, $dir );
    return $files;
}
