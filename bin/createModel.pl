#!/usr/bin/perl

=head1 NAME

createModel.pl

=head1 DESCRIPTION

This script creates a background model for a genome, computing the composition 
of the non-functional sequences and the elements presented naturally.

=head1 USAGE

perl createModel.pl -m MODEL [OPTIONS]

OPTIONS
    Parameter     Description                       Value          Default
    -h --help     Print this screen
    -v --verbose  Activate verbose mode

    -m --model    Create model                      String
    -d --dir      Write model in directory          Directory      ./data/MODEL
    -k --kmer     Profile k-mer size                Integer        4
    -w --win      Profile window size               Integer        200
    
    -f --fasta    Fasta sequences                   FileName*
    -r --repeats  RepeatMasker output               FileName*
    -t --trf      TRF output                        FileName*
    -g --genes    Gene annotation                   FileName*

  * File can be compressed (.gz/.bz2), if you are passing more than one file
    separate them with ',' Example: -f chr1.fa,chr2.fa,chr3.fa
    
=head1 EXAMPLES

You can obtain the models from our website, so you can we are currently suporting some
organisms with complete annotation in the UCSC Genome Database:
[http://hgdownload.cse.ucsc.edu/downloads.html]

You can create your own model fetching the data from the UCSC site:

  perl createModel.pl -m hg19

The model names are according to UCSC, like:
  - hg19    -> human genome release 19 (GRCh37)
  - hg18    -> human genome release 18 (GRCh36)
  - mm9     -> mouse genome release 9
  - equCab2 -> horse genome release 2

Th script will download the required files and process it. 

* Caution, the files are really big and the download time could be large.

Also you can use your own sequences and annotations to create a model:

  perl createModel.pl -m myOrg -f myOrg.fa -r myOrg_RM.out -t myOrg_TRF.out -g myOrg_Genes.table

In this case you must provide the sequences in a Fasta file, the RepeatMasker
output, the Tandem Repeat Finder output and the annotated genes in a tabular
text file with column 3 indicating the chromosome, column 5 the starting 
coordinate and column 6 indicating the ending coordinate.
  
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
my $model      = undef;      # Model definition
my $dir        = 'data';     # Output directory
my $kmer       = 4;          # K-mer size
my $win        = 200;        # Window size (non-overlapping)
my $fasta      = undef;      # Sequences in fasta files
my $repeat     = undef;      # RepeatMasker output
my $trf        = undef;      # TRF output
my $gene       = undef;      # Gene annotation
my $help       = undef;      # Help flag
my $verbose    = undef;      # Verbose mode flag
my $rm_tmp     = undef;      # Remove downloaded files

GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'm|model=s'        => \$model,
    'd|dir:s'          => \$dir,
    'k|kmer:i'         => \$kmer,
    'w|win:i'          => \$win,
    'f|fasta:s'        => \$fasta,
    'r|repeat:s'       => \$repeat,
    't|trf:s'          => \$trf,
    'g|gene:s'         => \$gene,
    'rm_tmp'           => \$rm_tmp
) or pod2usage(-verbose => 2);
pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $model);

# Configurable parameters
my $get         = 'wget -q'; # use this command to fetch files from internet
my $ucsc        = 'http://hgdownload.cse.ucsc.edu/goldenPath';
my $ucsc_genome = "$ucsc/$model/bigZips/chromFa.tar.gz";
my $ucsc_repeat = "$ucsc/$model/bigZips/chromOut.tar.gz";
my $ucsc_trf    = "$ucsc/$model/bigZips/chromTrf.tar.gz";
my $ucsc_gene   = "$ucsc/$model/database/ensGene.txt.gz"; 
