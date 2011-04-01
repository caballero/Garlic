#!/usr/bin/perl

=head1 NAME

createModel.pl

=head1 DESCRIPTION

This script creates a background model for a genome, computing the composition
of the non-functional sequences and the elements presented naturally.

=head1 USAGE

perl createModel.pl -m MODEL [OPTIONS]

OPTIONS
    Parameter          Description                   Value          Default
    -h --help          Print this screen
    -v --verbose       Activate verbose mode

    -m --model         Create model                  String
    -d --dir           Write model in directory      Directory      ./data
    -k --kmer          Profile k-mer size            Integer        4
    -w --win           Profile window size           Integer        200
    
    -f --fasta         Fasta sequences               FileName*
    -r --repeats       RepeatMasker output           FileName*
    -t --trf           TRF output                    FileName*
    -g --genes         Gene annotation               FileName*
    -e --exclude       Exclude this sequences        Pattern**
    
    --mask_repeat      Do intersersed repeats masking*** 
    --mask_trf         Do simple repeats masking***
    --no_mask_gene     Don't mask genes
    --no_repeat_table  Don't create repeats table
    --no_kmer_table    Don't compute k-mer composition

  * File can be compressed (.gz/.bz2), if you are passing more than one file 
    separate them with ',' Example: "-f chr1.fa,chr2.fa,chr3.fa"
 ** This pattern will match the fasta files. Example: "-e hap1" will exclude 
    all files with "hap1" in the name.
*** UCSC Genome Database sequences are repeat soft-masked.
    
=head1 EXAMPLES

You can obtain the models from our website, we are currently suporting some 
organisms with complete annotation in the UCSC Genome Database:
<http://hgdownload.cse.ucsc.edu/downloads.html>

You can create your own model fetching the data from the UCSC site:

  perl createModel.pl -m hg19

The model names are according to UCSC, like:
  - hg19    -> human genome release 19 (GRCh37)
  - hg18    -> human genome release 18 (GRCh36)
  - mm9     -> mouse genome release 9
  - equCab2 -> horse genome release 2

This script will download the required files and process it. Caution, the files
are really big and the download time could be large.

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

This is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation, 
either version 3 of the License, or (at your option) any later version.

This is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
code.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Find;

# Parameters initialization
my $model           = undef;      # Model definition
my $dir             = 'data';     # Output directory
my $kmer            = 4;          # K-mer size
my $win             = 200;        # Window size (non-overlapping)
my $fasta           = undef;      # Sequences in fasta files
my $repeat          = undef;      # RepeatMasker output
my $trf             = undef;      # TRF output
my $gene            = undef;      # Gene annotation
my $help            = undef;      # Help flag
my $verbose         = undef;      # Verbose mode flag
my $rm_tmp          = undef;      # Remove downloaded files
my $exclude         = undef;      # Exclude sequence name pattern
my $mask_repeat     = undef;      # Repeat masking flag
my $mask_trf        = undef;      # TRF masking flag
my $no_mask_gene    = undef;      # Gene masking flag
my $no_repeat_table = undef;      # No repeats profile flag
my $no_kmer_table   = undef;      # No kmer profile flag

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'm|model=s'       => \$model,
    'd|dir:s'         => \$dir,
    'k|kmer:i'        => \$kmer,
    'w|win:i'         => \$win,
    'f|fasta:s'       => \$fasta,
    'r|repeat:s'      => \$repeat,
    't|trf:s'         => \$trf,
    'g|gene:s'        => \$gene,
    'e|exclude:s'     => \$exclude,
    'rm_tmp'          => \$rm_tmp,
    'mask_repeat'     => \$mask_repeat,
    'mask_trf'        => \$mask_trf,
    'no_mask_gene'    => \$no_mask_gene,
    'no_repeat_table' => \$no_repeat_table,
    'no_kmer_table'   => \$no_kmer_table 
) or pod2usage(-verbose => 2);

# Call help if required
pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $model);

# Configurable parameters
my $get           = 'wget -c'; # command to fetch files from internet
   $get          .= ' -q' unless (defined $verbose);
my $unpack        = 'tar zxf'; # command to unpack the files downloaded
my $ucsc          = 'http://hgdownload.cse.ucsc.edu/goldenPath'; # UCSC url
my $ucsc_genome   = "$ucsc/$model/bigZips/chromFa.tar.gz";
my $ucsc_repeat   = "$ucsc/$model/bigZips/chromOut.tar.gz";
my $ucsc_trf      = "$ucsc/$model/bigZips/chromTrf.tar.gz";
my $ucsc_gene     = "$ucsc/$model/database/ensGene.txt.gz"; 

# Main variables
my @fasta    = ();
my @repeat   = ();
my @trf      = ();
my @gene     = ();
my %seq      = ();
my %data     = ();
my %kmer     = ();
my ($seq, $seq_id, $ini, $end, $len); 
my @gc = qw/0-5 5-10 10-15 15-20 20-25 25-30 30-35 35-40 40-45 45-50 50-55 55-60 65-70 70-75 75-80 80-85 85-90 90-95 95-100/;

# Check directories, create them if required
unless (-e "$dir" and -d "$dir") {
    warn "creating directory $dir\n" if (defined $verbose);
    mkdir "$dir";
}
unless (-e "dir/$model" and -d "$dir/$model") {
    warn "creating directory $dir/$model\n" if (defined $verbose);
    mkdir "$dir/$model";
}
warn "moving to $dir/$model\n" if (defined $verbose);
chdir "$dir/$model" or die "cannot move to $dir/$model";

# Get files from UCSC genome database if required
getUCSC_fasta()  unless (defined $fasta);
getUCSC_repeat() unless (defined $repeat);
getUCSC_trf()    unless (defined $trf);
getUCSC_gene()   unless (defined $gene);

# Creating list of files to process
@fasta   = split (/,/, $fasta);
@repeat  = split (/,/, $repeat);
@trf     = split (/,/, $trf);
@gene    = split (/,/, $gene);
$exclude =~ s/,/|/g;

# Load sequences and mask them
readFasta();
maskRepeat()   if (defined $mask_repeat);
maskTRF()      if (defined $mask_trf);
maskGene() unless (defined $no_mask_gene);

# Create the K-mer/Window table
profileSeq() unless (defined $no_kmer_table);

# Create the Repeat Table
profileRepeat() unless (defined $no_repeat_table);

removeTmp() if (defined $rm_tmp);

warn "Done\n" if (defined $verbose);


#################################################
##      S  U  B  R  O  U  T  I  N  E  S        ##
#################################################

sub searchFiles {
    my ($pat, $dir) = @_;
    my $files = undef;
    find ( sub { $files .= $File::Find::name . ',' if (m/$pat/) }, $dir );
    return $files;
}

sub getUCSC_gene {
    warn "obtaining Gene files from $ucsc\n" if (defined $verbose);
    system ("$get $ucsc_gene");
    die "cannot find Gene output in $ucsc_gene" unless (-e 'ensGene.txt.gz' and -s 'ensGene.txt.gz');
    $gene = 'ensGene.txt.gz';
}

sub getUCSC_trf {
    warn "obtaining TRF files from $ucsc\n" if (defined $verbose);
    mkdir 'TRF' unless (-e 'TRF' and -d 'TRF');
    chdir 'TRF' or die "cannot move to TRF directory\n";
    system ("$get $ucsc_trf");
    die "cannot find TRF output in $ucsc_trf" unless (-e 'chromTrf.tar.gz' and -s 'chromTrf.tar.gz');
    warn "unpacking TAR\n" if (defined $verbose);
    system ("$unpack chromTrf.tar.gz");
    chdir '..';
    warn "searching TRF files\n" if (defined $verbose);
    $trf = searchFiles('.bed$', 'TRF');
}

sub getUCSC_repeat {
    warn "obtaining RepeatMasker files from $ucsc\n" if (defined $verbose);
    mkdir 'RM' unless (-e 'RM' and -d 'RM');
    chdir 'RM' or die "cannot move to RM directory\n";;
    system ("$get $ucsc_repeat");
    die "cannot find RepeatMasker output in $ucsc_repeat" unless (-e 'chromOut.tar.gz' and -s 'chromOut.tar.gz');
    warn "unpacking TAR\n" if (defined $verbose);
    system ("$unpack chromOut.tar.gz");
    chdir '..';
    warn "searching RM files\n" if (defined $verbose);
    $repeat = searchFiles('.out$', 'RM');
}

sub getUCSC_fasta {
    warn "obtaining fasta files from $ucsc\n" if (defined $verbose);
    mkdir 'fasta' unless (-e 'fasta' and -d 'fasta');
    chdir 'fasta' or die "cannot move to fasta directory\n";
    system ("$get $ucsc_genome");
    die "cannot find genomic sequences in $ucsc_genome" unless (-e 'chromFa.tar.gz' and -s 'chromFa.tar.gz');
    warn "unpacking TAR\n" if (defined $verbose);
    system ("$unpack chromFa.tar.gz");
    chdir '..';
    warn "searching fasta files\n" if (defined $verbose);
    $fasta = searchFiles('.fa$', 'fasta');
}

sub readFasta {
    warn "loading fasta sequences\n" if (defined $verbose);
    foreach my $file (@fasta) {
        if (defined $exclude) {
            next if ($file =~ m/$exclude/);
        }
        warn "  reading $file\n" if (defined $verbose);
        my $fileh = $file;
        $fileh = "gunzip  -c $file | " if ($file =~ m/\.gz$/);
        $fileh = "bunzip2 -c $file | " if ($file =~ m/\.bz2$/);
        open FH, "$fileh" or die "cannot open $file\n";
        while (<FH>) {
            chomp;
            if (m/^>(\S+?)/) { $seq_id = $1;        }
            else             { $seq{$seq_id} .= $_; }
        }
        close FH;
    }
}

sub maskRepeat {
    warn "masking repeats\n" if (defined $verbose);
    foreach my $file (@repeat) {
        my $fileh = $file;
        $fileh = "gunzip  -c $file | " if ($file =~ m/\.gz$/);
        $fileh = "bunzip2 -c $file | " if ($file =~ m/\.bz2$/);
        open FH, "$fileh" or die "cannot open $file\n";
        while (<FH>) {
            s/^\s*//;
            my @arr = split (/\s+/, $_);
            $seq_id = $arr[4];
            next unless (defined $seq{$seq_id});
            $ini = $arr[5];
            $end = $arr[6];
            $len = $end - $ini;
            $ss  = substr ($seq{$seq_id}, $ini - 1, $len);
            substr ($seq{$seq_id}, $ini - 1, $len) = lc $ss;
        }
        close FH;
    }
}

sub maskTRF {
    warn "masking simple repeats\n" if (defined $verbose);
    foreach my $file (@trf) {
        my $fileh = $file;
        $fileh = "gunzip  -c $file | " if ($file =~ m/\.gz$/);
        $fileh = "bunzip2 -c $file | " if ($file =~ m/\.bz2$/);
        open FH, "$fileh" or die "cannot open $file\n";
        while (<FH>) {
            my @arr = split (/\t/, $_);
            $seq_id = $arr[0];
            next unless (defined $seq{$seq_id});
            $ini = $arr[1];
            $end = $arr[2];
            $len = $end - $ini;
            $ss  = substr ($seq{$seq_id}, $ini - 1, $len);
            substr ($seq{$seq_id}, $ini - 1, $len) = lc $ss;
        }
        close FH;
    }
}

sub maskGenes {
    warn "masking genes\n" if (defined $verbose);
    foreach my $file (@gene) {
        my $fileh = $file;
        $fileh = "gunzip  -c $file | " if ($file =~ m/\.gz$/);
        $fileh = "bunzip2 -c $file | " if ($file =~ m/\.bz2$/);
        open FH, "$fileh" or die "cannot open $file\n";
        while (<FH>) {
            my @arr = split (/\t/, $_);
            $seq_id = $arr[2];
            next unless (defined $seq{$seq_id});
            $ini = $arr[4];
            $end = $arr[5];
            $len = $end - $ini;
            $ss  = substr ($seq{$seq_id}, $ini - 1, $len);
            substr ($seq{$seq_id}, $ini - 1, $len) = lc $ss;
        }
        close FH;
    }
}

sub profileSeq {
    warn "profiling sequences, k-mer=$kmer and window=$win\n" if (defined $verbose);
    my $bp_slices = 0;
    createKmer($kmer - 1);
    foreach $seq_id (keys %seq) {
	warn "  analyzing sequence $seq_id\n" if (defined $verbose);
        $seq = $seq{$seq_id};
	$seq =~ s/[^ACGT]/N/g;
	$seq =~ s/N+/N/g;
	for (my $i = 0; $i <= (length $seq) - $win; $i++) {
	    my $slice = substr ($seq, $i, $win);
	    next if ($slice =~ m/N/);
	    $bp_slices += $win;
	    my $gc = calcGC($slice);
	    for (my $j = 0; $j <= $win - $kmer; $j++) {
	        my $kmer = substr ($slice, $j, $kmer);
		$data{$gc}{$kmer}++;
	    }
        }
    }
    warn "  $bp_slices bases analyzed\n" if (defined $verbose);
    my $file = "$model.K$kmer.W$win.data";
    warn "writing k-mer profile in $file\n" if (defined $verbose);
    open K, ">$file" or die "cannot write file $file\n";
    foreach my $gc (@gc) {
	print K "#GC=$gc\n";
	foreach my $word (@kmer) {
            my $tot = 0;
	    foreach my $b (@bases) {
                $tot += $data{$gc}{"$word$b"} if (defined $data{$gc}{"$word$b"});
	    }
	    if ($tot > 0) {
		foreach my $b (@bases) {
		    my $cnt = 0;
		    $cnt = $data{$gc}{"$word$b"} if (defined $data{$gc}{"$word$b"});
		    $frq = sprintf ("%.8f", $cnt / $tot);
		    print K "$word$b\t$frq\t$cnt\n";
	        }
	    }
	    else {
		foreach my $b (@bases) {
		    print K "$word$b\t0.25\t0\n";
	        }
	    }
        }
    }
}

sub createKmer {

}

sub profileRepeat{
    warn "profiling repeats\n" if (defined $verbose);
    
}

sub removeTmp {
    warn "removing temporary files\n" if (defined $verbose);
    my @files = @gene;
    push @files, 'fasta';
    push @files, 'RM';
    push @files, 'TRF';
    foreach my $file (@files) {
        system ("rm -rf $file"); 
    }
}

sub calcGC {
    my $seq = shift @_;
    my $len = $seq =~ tr/ACGTacgt/ACGTacgt/;
    my $ngc = $seq =~ tr/CGcg/CGcg/;
    return 'NA' if ($len < 1);

    my $gc  = $ngc / $len;
    if    ($gc <  5) { $gc =  '0-5'  ; }
    elsif ($gc < 10) { $gc =  '5-10' ; }
    elsif ($gc < 15) { $gc = '10-15' ; }    
    elsif ($gc < 20) { $gc = '15-20' ; }   
    elsif ($gc < 25) { $gc = '20-25' ; }
    elsif ($gc < 30) { $gc = '25-30' ; }
    elsif ($gc < 35) { $gc = '30-35' ; }
    elsif ($gc < 40) { $gc = '35-40' ; }
    elsif ($gc < 45) { $gc = '40-45' ; }
    elsif ($gc < 50) { $gc = '45-50' ; }
    elsif ($gc < 55) { $gc = '50-55' ; }
    elsif ($gc < 60) { $gc = '55-60' ; }   
    elsif ($gc < 65) { $gc = '60-65' ; }
    elsif ($gc < 70) { $gc = '65-70' ; }
    elsif ($gc < 75) { $gc = '70-75' ; }
    elsif ($gc < 80) { $gc = '75-80' ; }
    elsif ($gc < 85) { $gc = '80-85' ; }
    elsif ($gc < 90) { $gc = '85-90' ; }
    elsif ($gc < 95) { $gc = '90-95' ; }
    else             { $gc = '95-100'; }
    
    return $gc;
}
    
