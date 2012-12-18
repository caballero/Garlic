#!/usr/bin/perl

=head1 NAME

createSequence_web.pl

=head1 DESCRIPTION

Web form for createFakeSequence.pl

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
use CGI::Pretty qw/:standard/;

my $create   = './createFakeSequence.pl';
my $compress = 'zip';
my $min_size = 1000;
my $max_size = 100000;
my @gc       = (0,10,20,30,40,50,60,70,80,90,100);
my @models   = qw/hg19 ailMel1 bosTau4 equCab2 gasAcu1 mm9 oryLat2 susScr2 anoCar2 calJac3 ci2 felCat4 monDom5 oviAri1 rheMac2 canFam2 danRer7 fr2 loxAfr3 ornAna1 panTro3 rn4 tetNig2 apiMel2 cavPor3 dm3 galGal3 mm10 oryCun2 strPur2	xenTro2/;

# CSS definition
my $style    =<<__STYLE__
<style type="text/css">
    body  { color: navy; background: lightyellow }
    h1    { color: navy; }
    table { border: 2px solid grey; border-collapse: collapse }
    th    { text-align: center; font-weight: bold; border: 2px solid grey; padding: 5px }
    td    { text-align: left; border: 1px solid grey; padding: 5px }
</style>
__STYLE__
;

# HTML header
print header('text/html'); 
print "<html>\n<head>\n<title></title>\n$style\n</head>\n<body>";
print h2('Artificial Intergenic Sequence Generator'), hr();
print p('This is a simulator designed to produce a DNA sequence similar to the non-functional (intergenic) regions of a genome to be use in genomic analysis.');
print p('If you are planning to run a large scale simulations, please install a local version, the source code is available in <a href="http://caballero.github.com/FakeSequence/">here </a>and the models <a href="models/">here.</a>');
print hr();

if (defined param('model')) {
    my $model   = param('model');
    my $size    = param('size');
    my $kmer    = param('kmer');
    my $win     = param('win');
    my $min_gc  = param('min_gc');
    my $max_gc  = param('max_gc');
    my $repfrac = param('repfrac');
    my $reptype = param('reptype');
    
    if ($size < $min_size or $size > $max_size) {
        print p(i("Please use a size between $min_size and $max_size bases"));
        createForm();
    }
    else {
        print p('Generating sequences, please wait ...');
        my $id  = int(rand(1e9));
        my $cmd = "$create -m $model -s $size -n tmp/$id -k $kmer -r $repfrac -g $min_gc -c $max_gc ";
        $cmd   .= "-t $reptype" if ($reptype ne 'All');
        system ($cmd);
        system ("$compress tmp/$id.zip tmp/$id.*");
        system ("rm tmp/$id.fasta tmp/$id.inserts");
        print p("Your files are in <a href\"tmp/$id.zip\">$id.zip</a>"); 
    }
}
else {
    createForm();
}

# HTML footer
print p('Citation: Caballero  ...');
print p('Contact: please use the <a href="https://github.com/caballero/FakeSequence/issues">GitHub form.</a>');
print p(i('Institute for Systems Biology (2012)'));
print end_html();

sub createForm {
    print start_form();
    print p('Model ', popup_menu(-name => 'model', -values => \@models, -default => 'hg19'));
    print p('Kmer ',  popup_menu(-name => 'kmer', -values => [4, 6, 8], -default => 4));
    print p('Window ',  popup_menu(-name => 'win', -values => [1000], -default => 1000));
    print p('Size ',  textfield(-name => 'size', -value => 10000, -size => 10), i(" number between $min_size - $max_size"));
    print p('Minimal G+C ', popup_menu(-name => 'min_gc', -values => \@gc, -default => 40));
    print p('Maximal G+C ', popup_menu(-name => 'max_gc', -values => \@gc, -default => 60));
    print p('Repetitive fraction ',  textfield(-name => 'repfrac', -value => 0.5, -size => 10), i(" number between 0.0 - 0.95"));
    print p('Repeat Type ', textfield(-name => 'reptype', -value => 'All', -size => 10), i(" like Alu, MIR, LINE, LTR, ..."));

    print br();
	print submit(-name => 'Create!');
	print end_form();
}

