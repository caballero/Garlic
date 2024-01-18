#!/usr/local/bin/perl

=head1 NAME

createModel.pl

=head1 DESCRIPTION

This script creates a background model for a genome, computing the composition
of the non-functional, non-repetitive sequences and the repetitive elements 
presented naturally.

=head1 USAGE

perl createModel.pl -m MODEL [OPTIONS]

OPTIONS
    Parameter          Description                   Value          Default
    -h --help          Print this screen
    -v --verbose       Activate verbose mode

    -m --model         Create model                  String
    -d --dir           Write model in directory      Directory      ./data
    -k --kmer          Profile k-mer size            Integer        4
    -w --win           Profile window size           Integer        1000
    -b --binsize       Bin size (for %GC)            Integer        10000
    -n --numgc         Number of GC bins             Integer        25
    
    -f --fasta         Fasta sequences               FileName*
    -r --repeats       RepeatMasker output           FileName*
    -t --trf           TRF output                    FileName*
    -g --genes         Gene annotation               FileName*
    -a --analyze       Analyze these regions         FileName*
    -e --exclude       Exclude these sequences       Pattern**
    
    --mask_exon        Mask exons (default is whole genes)
    --no_intron_filter Don't filter intronic sequences
    --no_mask_gene     Don't mask genes
    --no_mask_repeat   Don't mask repeats
    --no_mask_trf      Don't mask simple repeats
    --no_repeat_table  Don't create repeats table
    --no_kmer_table    Don't compute k-mer composition
    --rm_tmp           Remove temp files
    --keep_dw_files    Keep downloaded files (overrules --rm_tmp)
    --gc_post_mask     Compute GC bins after masking (default is before)
    --revcomp_kmer     Count kmers in reverse complement chain
    --write_mask_seq   Write the masked sequence

  * File can be compressed (.gz/.bz2), if you are passing more than one file 
    separate them with commas. Example: "-f chr1.fa,chr2.fa,chr3.fa"
 ** This pattern will match the fasta files. Example: "-e hap" will exclude 
    all files with "hap" in the name.
    
=head1 EXAMPLES

You can obtain the models from our website, we are currently suporting some 
organisms with complete annotation in the UCSC Genome Database:
<http://hgdownload.cse.ucsc.edu/downloads.html>

You can create your own model fetching the data from the UCSC site:

  perl createModel.pl -m hg19

The model names are according to UCSC nomenclature:
  - hg19    -> human genome release 19 (GRCh37)
  - hg18    -> human genome release 18 (GRCh36)
  - mm9     -> mouse genome release 9
  - equCab2 -> horse genome release 2
  - ...

This script will download the required files and process it. Caution, some files
are really big and the download time may be large.

Alternatively, you can use your own sequences and annotations to create a model:

  perl createModel.pl -m myOrg -f myOrg.fa -r RM.out -t TRF.out -g Genes.table

In this case you must provide the sequences in a Fasta file, the RepeatMasker 
output, the Tandem Repeat Finder output and the annotated genes in a tabular 
text file similar to the UCSC ensGenes.txt format (not knownGene.txt).
  
=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

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
use Cwd;

# Parameters initialization
my $model           = undef;     # Model definition
my $dir             = 'data';    # Output directory
my $kmer            = 4;         # K-mer size
my $win             = 1000;      # Window size (non-overlapping)
my $binsize         = 10000;     # Bin size
my $gcBuckets       = 25;        # Number of GC buckets/bins
my $fasta           = undef;     # Sequences in fasta files
my $repeat          = undef;     # RepeatMasker output
my $trf             = undef;     # TRF output
my $gene            = undef;     # Gene annotation
my $region          = undef;     # Define a region to analize
my $help            = undef;     # Help flag
my $verbose         = undef;     # Verbose mode flag
my $rm_tmp          = undef;     # Remove temporal files
my $exclude         = undef;     # Exclude sequences by name pattern
my $no_mask_gene    = undef;     # Gene masking flag
my $no_mask_repeat  = undef;     # Repeats masking flag
my $no_mask_trf     = undef;     # TRF masking flag
my $no_repeat_table = undef;     # No repeats profile flag
my $no_kmer_table   = undef;     # No kmer profile flag
my $write_mask_seq  = undef;     # Write fasta flag
my $keep_dw_files   = undef;     # Keep downloaded files
my $gc_post_mask    = undef;     # Compute GC bins before masking
my $revcomp_kmer    = undef;     # Count kmers in the reverse-complement chain
my $mask_exon       = undef;     # Mask exonic regions (default is whole gene)
my $no_intron       = undef;     # No intron filter flag

# Fetch options
GetOptions(
            'h|help'           => \$help,
            'v|verbose'        => \$verbose,
            'm|model=s'        => \$model,
            'n|numgc:i'        => \$gcBuckets,
            'd|dir:s'          => \$dir,
            'k|kmer:i'         => \$kmer,
            'w|win:i'          => \$win,
            'b|binsize:i'      => \$binsize,
            'f|fasta:s'        => \$fasta,
            'r|repeat:s'       => \$repeat,
            't|trf:s'          => \$trf,
            'g|gene:s'         => \$gene,
            'e|exclude:s'      => \$exclude,
            'a|analyze:s'      => \$region,
            'rm_tmp'           => \$rm_tmp,
            'no_mask_gene'     => \$no_mask_gene,
            'no_mask_repeat'   => \$no_mask_repeat,
            'no_mask_trf'      => \$no_mask_trf,
            'no_repeat_table'  => \$no_repeat_table,
            'no_kmer_table'    => \$no_kmer_table,
            'write_mask_seq'   => \$write_mask_seq,
            'keep_dw_files'    => \$keep_dw_files,
            'gc_post_mask'     => \$gc_post_mask,
            'revcomp_kmer'     => \$revcomp_kmer,
            'mask_exon'        => \$mask_exon,
            'no_intron_filter' => \$no_intron
    )
    or pod2usage( -verbose => 2 );

# Call help if required
pod2usage( -verbose => 2 ) if ( defined $help );
pod2usage( -verbose => 2 ) if !( defined $model );

# UCSC model files
my %files = ();
loadFiles();

# Configurable parameters
my $get = 'wget -c';    # command to fetch files from internet
$get .= ' -q' unless ( defined $verbose );
my $unpack = 'tar zxf';      # command to unpack the files downloaded
my $unzip  = 'unzip';        # command to decompress the files downloaded
my $gunzip = 'gunzip -c';    # command to decompress the files downloaded
my $ucsc        = 'http://hgdownload.cse.ucsc.edu/goldenPath';        # UCSC url
my $repeatMasker= 'http://www.repeatmasker.org/genomes';               # RepeatMasker website url
#my $ucsc_repeat = "$ucsc/$model/bigZips/" . $files{$model}{'RMA'};

## TODO: This should be deprecated.  The main program should always expect the user
##       to provide these files as parameters.  A seperate utility could be created
##       to download these files from various sources if desired, although maintaining
##       that in the long run may not be worthwhile.
my $ucsc_genome;
my $rm_repeat;
my $ucsc_trf;
my $ucsc_gene;
if ( exists $files{$model} ) {
  $ucsc_genome = "$ucsc/$model/bigZips/" . $files{$model}{'FAS'};
  $rm_repeat   = "$repeatMasker/" . $files{$model}{'RMA'};
  $ucsc_trf    = "$ucsc/$model/bigZips/" . $files{$model}{'TRF'};
  $ucsc_gene   = "$ucsc/$model/database/" . $files{$model}{'GEN'};
}

# Main variables
my @fasta     = ();            # list of fasta files
my @repeat    = ();            # list of RepeatMasker output files
my @trf       = ();            # list of TRF ouput files
my @gene      = ();            # list of gene annotation files
my @region    = ();            # list of selected regions
my %bingc     = ();            # index to search the regional GC
my %seq       = ();            # complete sequences
my %repeat    = ();            # repeats information
my %repinsert = ();            # repeat insertions
my %suminsert = ();            # total sum in repeat insertions
my %genes     = ();            # gene coordinates
my %region    = ();            # preselected regions
my %ebases    = ();            # effective bases per GC bin
my %rbases    = ();            # interspersed repeats bases per GC bin
my %sbases    = ();            # low complexity/simple repeat bases per GC bin
my @dna       = qw/A C G T/;
my ( $seq, $ss, $seq_id, $ini, $end, $len, $name );
my @gc     = ();               # GC content bins
my @gcBins = ();               # cutoffs for content bins

# Check directories, create them if required
unless ( -e "$dir" and -d "$dir" )
{
  warn "creating directory $dir\n" if ( defined $verbose );
  mkdir "$dir";
}
unless ( -e "$dir/$model" and -d "$dir/$model" )
{
  warn "creating directory $dir/$model\n" if ( defined $verbose );
  mkdir "$dir/$model";
}
warn "moving to $dir/$model\n" if ( defined $verbose );

my $startDir = cwd();

# Get files from UCSC genome database if required
print "Getting data from UCSC if required...\n";

my $outDir = "$dir/$model";

sub buildAbsPathFileList {
  my $fileStr = shift;

  my @files;
  foreach my $file ( split(/,/, $fileStr) ) {
    # Convert relative paths (if used) to absolute
    push @files, Cwd::abs_path($file);
  }
  return @files;
}

if ( defined $fasta && $fasta ne "" ) {
  @fasta = buildAbsPathFileList($fasta);
}else {
  getUCSC_fasta($outDir);
}

unless( $no_repeat_table ) {
  if ( defined $repeat && $repeat ne "" ) {
    @repeat = buildAbsPathFileList($repeat);
  }else {
    getUCSC_repeat($outDir);
  }
  if ( defined $trf && $trf ne "" ) {
    @trf = buildAbsPathFileList($trf);
  }else {
    getUCSC_trf($outDir);
  }
}

if ( defined $gene && $gene ne "" ) {
  @gene = buildAbsPathFileList($gene);
}else {
  getUCSC_gene($outDir);
}

if ( defined $region && $region ne "" ) {
  @region = buildAbsPathFileList($region);
}

$exclude =~ s/,/|/g if ( defined $exclude );

# The remainder of the program assumes we are working in $dir/$model
chdir "$dir/$model" or die "cannot move to $dir/$model";

# Load sequences and mask it
print "Loading sequences...\n";
readFasta();
print "Calc GC bins...\n";
calcBinGC() unless ( defined $gc_post_mask );

unless ( defined $no_mask_repeat ) {
  print "Masking repeats...\n";
  maskRepeat();
}

unless ( defined $no_mask_trf ) {
  print "Masking TRF...\n";
  maskTRF();
}

if ( defined $mask_exon )
{
  print "Masking Exons...\n";
  maskExon();
} else
{
  unless ( defined $no_mask_gene ) {
    print "Masking Genes...\n";
   maskGene();
  }
}

calcBinGC() if ( defined $gc_post_mask );

if ( defined $write_mask_seq )
{
  print "Writing masked seq...\n";
  writeMaskSeq( "$model.masked.fa" );
}

# Create the K-mer/Window and GCt tables
unless ( defined $no_kmer_table )
{
  print "Creating K-mer/Window and GCt tables...\n";
  profileSeqs();
}

# Create the Repeat Table
print "Creating the repeats table\n";
profileRepeats() unless ( defined $no_repeat_table );

# Clean up
removeTmp() if ( defined $rm_tmp );

warn "Done\n" if ( defined $verbose );

#################################################
##      S  U  B  R  O  U  T  I  N  E  S        ##
#################################################

sub searchFiles
{

  # deep search of files/patterns in a directory
  my ( $pat, $dir ) = @_;
  my $files = undef;
  find( sub { $files .= $File::Find::name . ',' if ( m/$pat/ ) }, $dir );
  return $files;
}

sub getUCSC_gene
{
  my $outDir = shift;

  # Gene models to mask functional regions
  $gene = $files{$model}{'GEN'};
  die "Sorry, $model gene annotation is missing\n" unless ( defined $gene );

  if ( -e $gene )
  {
    warn "File $gene present, using it\n" if ( defined $verbose );
  } else
  {
    warn "obtaining Gene files from $ucsc\n" if ( defined $verbose );
    system( "$get --directory-prefix=$outDir $ucsc_gene" );
    die "cannot find gene ann in $ucsc_gene\n" unless ( -e $gene );
  }
}

sub getUCSC_trf
{
  my $outDir = shift;

  # TRF output is used for simple repeat annotation
  my $target_file = $files{$model}{'TRF'};
  die "Sorry, $model TRF file is missing\n" unless ( defined $target_file );

  if ( -e "TRF/$target_file" )
  {
    warn "$target_file present, using it\n" if ( defined $verbose );
    chdir 'TRF' or die "cannot move to TRF directory\n";

  } else
  {
    warn "obtaining TRF files from $ucsc\n" if ( defined $verbose );
    mkdir 'TRF' unless ( -e 'TRF' and -d 'TRF' );
    chdir 'TRF' or die "cannot move to TRF directory\n";
    system( "$get --directory-prefix=$outDir $ucsc_trf" );
    die "cannot find TRF output in $ucsc_trf\n" unless ( -e $target_file );
  }

  unpackFiles( $target_file );

  chdir '..';
  warn "searching TRF files\n" if ( defined $verbose );
  $trf = searchFiles( '.bed$', 'TRF' );
  die "cannot find TRF files!" unless ( defined $trf );
}

sub getUCSC_repeat
{
  my $outDir = shift;
  # RepeatMasker output
  my $target_file = $files{$model}{'RMA'};

  if ( $target_file =~ /^\s*\S+\/\S+\/(\S+\.fa.align.gz)\s*$/ )
  {
    $target_file = $1;
  }
  die "Sorry, $model RepeatMasker file is missing\n"
      unless ( defined $target_file );

  if ( -e "RM/$target_file" )
  {
    warn "$target_file present, using it\n" if ( defined $verbose );
    chdir 'RM' or die "cannot move to RM directory\n";

  } else
  {
    #warn "obtaining RepeatMasker files from $ucsc\n" if ( defined $verbose );
    warn "obtaining RepeatMasker files from $repeatMasker\n" if ( defined $verbose );
    mkdir 'RM' unless ( -e 'RM' and -d 'RM' );
    chdir 'RM' or die "cannot move to RM directory\n";
    system( "$get --directory-prefix=$outDir $rm_repeat" );
    die "cannot find RepeatMasker align in $rm_repeat\n"
        unless ( -e $target_file );
  }

  unpackFiles( $target_file );

  chdir '..';
  warn "searching RM files\n" if ( defined $verbose );
  $repeat = searchFiles( '.align$', 'RM' );
  die "cannot find RM files!" unless ( defined $repeat );
}

sub getUCSC_fasta
{
  my $outDir = shift;

  # Chromosomal sequences
  my $target_file = $files{$model}{'FAS'};
  die "Sorry, $model don't have fasta\n" unless ( defined $target_file );

  if ( -e "fasta/$target_file" )
  {
    warn "$target_file present, using it\n" if ( defined $verbose );
    chdir 'fasta' or die "cannot move to fasta directory\n";
  } else
  {
    warn "obtaining fasta files from $ucsc\n" if ( defined $verbose );
    mkdir 'fasta' unless ( -e 'fasta' and -d 'fasta' );
    chdir 'fasta' or die "cannot move to fasta directory\n";
    system( "$get --directory-prefix=$outDir $ucsc_genome" );
    die "cannot find genomic seq in $ucsc_genome\n" unless ( -e $target_file );
  }

  unpackFiles( $target_file );

  chdir '..';
  warn "searching fasta files\n" if ( defined $verbose );
  $fasta = searchFiles( '.fa$', 'fasta' );
  die "cannot find fasta files!" unless ( defined $fasta );
}

sub unpackFiles
{
  my $file = shift;
  if ( $file =~ m/tar.gz$/ )
  {
    warn "unpacking TAR: $file\n" if ( defined $verbose );
    system( "$unpack $file" );
  } elsif ( $file =~ m/.zip$/ )
  {
    warn "unpacking ZIP: $file\n" if ( defined $verbose );
    system( "$unzip $file" );
  } elsif ( $file =~ m/.gz$/ )
  {
    warn "unzipping GZ: $file\n" if ( defined $verbose );
    my $out_file = $file;
    $out_file =~ s/.gz$//;
    system( "$gunzip $file > $out_file" );
  } else
  {
    die "Abort, cannot recognize the file format of $file\n";
  }
}

sub defineFH
{

  # define how to read a file
  my $fi = shift;
  my $fh = $fi;
  $fh = "gunzip  -c $fi | " if ( $fi =~ m/\.gz$/ );
  $fh = "bunzip2 -c $fi | " if ( $fi =~ m/\.bz2$/ );
  return $fh;
}

sub readFasta
{

  # read tha fasta files, apply filters if required
  warn "loading fasta sequences\n" if ( defined $verbose );
  foreach my $file ( @fasta )
  {
    if ( defined $exclude )
    {
      next if ( $file =~ m/$exclude/ );
    }
    warn "    reading $file\n" if ( defined $verbose );
    my $fileh = defineFH( $file );
    open FH, "$fileh" or die "cannot open $fileh in " . cwd() . "\n";
    while ( <FH> )
    {
      chomp;
      if ( /^>(\S+)/ )
      {
        $seq_id = $1;
      } else
      {
        s/[^ACGTNacgt]/N/g;
        $seq{$seq_id} .= $_;
      }
    }
    close FH;
  }
}

# Should work with both *.align and *.out files
sub maskRepeat
{

  # mask sequences with RepeatMasker annotation
  warn "masking repeats\n" if ( defined $verbose );
  my %mask = ();
  foreach my $file ( @repeat )
  {
    if ( defined $exclude )
    {
      next if ( $file =~ m/$exclude/ );
    }
    my $fileh = defineFH( $file );
    open FH, "$fileh" or die "cannot open >$fileh<: $!\n";
    while ( <FH> )
    {
      s/^\s*//;
      next unless ( m/^\d+\s+\d+\.\d+/ );    # skip headers
      next if ( m/Simple_repeat|Low_complexity|Unknown|Satellite/ );
      my @arr = split( /\s+/, $_ );
      $seq_id = $arr[ 4 ];
      next unless ( exists $seq{$seq_id} );
      $ini = $arr[ 5 ];
      $end = $arr[ 6 ];
      $len = $end - $ini - 1;
      if ( $len > 0 ) {
        substr( $seq{$seq_id}, $ini - 1, $len ) = 'R' x $len;
      }
    }
    close FH;
  }
}

sub maskTRF
{

  # mask sequences with TRF annotation
  warn "masking simple repeats\n" if ( defined $verbose );
  my %mask = ();
  my $ranges = 0;
  foreach my $file ( @trf )
  {
    if ( defined $exclude )
    {
      next if ( $file =~ m/$exclude/ );
    }
    my $fileh = defineFH( $file );
    open FH, "$fileh" or die "cannot open $file\n";
    while ( <FH> )
    {
      my @arr = split( /\s+/, $_ );
      $seq_id = $arr[ 0 ];
      if ( ! exists $seq{$seq_id} )
      {
	warn "Simple repeat target sequence \'$seq_id\' is not in the genome!\n";
	next;
      }
      $ini = $arr[ 1 ];
      $end = $arr[ 2 ];
      $len = $end - $ini - 1;
      # RMH: This is an inconsistency.  TRF output is normally 1-based/fully-closed, however,
      #      the input is assuming TRF data from UCSC which is in bed format 0-based/half-open.
      #      The original code subtracted 1 from the start which will make the substr function
      #      wrap around to the end.  The new line assumes 0-based coordinates and leaves the
      #      start alone.
      #substr( $seq{$seq_id}, $ini - 1, $len ) = 'S' x $len;
      substr( $seq{$seq_id}, $ini, $len ) = 'S' x $len;
      $ranges++;
    }
    close FH;
  }
  if ( $ranges < 1 )
  {
    print "  - There were no simple repeats to mask.\n";
  }
  warn "  - Masked $ranges simple repeat ranges.\n" if ( defined $verbose);
}

sub maskExon
{

  # mask exon sequences with gene annotation
  warn "masking exons\n" if ( defined $verbose );
  my %mask = ();
  foreach my $file ( @gene )
  {
    my $fileh = defineFH( $file );
    open FH, "$fileh" or die "cannot open $file\n";
    while ( <FH> )
    {
      my @arr = split( /\s+/, $_ );
      $seq_id = $arr[ 2 ];
      next unless ( defined $seq{$seq_id} );
      $arr[ 9 ]  =~ s/,$//;
      $arr[ 10 ] =~ s/,$//;
      my @ini = split( /,/, $arr[ 9 ] );
      my @end = split( /,/, $arr[ 10 ] );
      for ( my $i = 0 ; $i <= $#ini ; $i++ )
      {
        $len = $end[ $i ] - $ini[ $i ] - 1;
        substr( $seq{$seq_id}, $ini[ $i ] - 1, $len ) = 'X' x $len;

      }
    }
    close FH;
  }
}

sub maskGene
{

  # mask sequences with gene annotation
  warn "masking genes\n" if ( defined $verbose );
  loadGenes();
  foreach $seq_id ( keys %genes )
  {
    next unless ( defined $seq{$seq_id} );
    foreach my $block ( @{ $genes{$seq_id} } )
    {
      ( $ini, $end, $name ) = split( /\s+/, $block );
      $len = $end - $ini - 1;
      substr( $seq{$seq_id}, $ini, $len ) = 'X' x $len;
    }
  }
  %genes = ();
}

sub loadGenes
{

  # load gene annotation
  warn "reading coordinates for genes\n" if ( defined $verbose );
  foreach my $file ( @gene )
  {
    my $fileh = defineFH( $file );
    open FH, "$fileh" or die "cannot open $file\n";
    while ( <FH> )
    {
      my @arr = split( /\s+/, $_ );
      $seq_id = $arr[ 2 ];
      next unless ( defined $seq{$seq_id} );
      $ini = $arr[ 4 ];
      $end = $arr[ 5 ];
      push @{ $genes{$seq_id} }, "$ini\t$end\tgene";
    }
    close FH;
  }
}

sub loadGenesBin
{

  # load gene annotation
  warn "reading coordinates for genes\n" if ( defined $verbose );
  my $bin_size = $binsize * 10;
  my $bin_     = undef;
  foreach my $file ( @gene )
  {
    my $fileh = defineFH( $file );
    open FH, "$fileh" or die "cannot open $file\n";
    while ( <FH> )
    {
      my @arr = split( /\s+/, $_ );
      $seq_id = $arr[ 2 ];
      next unless ( defined $seq{$seq_id} );
      $ini  = $arr[ 4 ];
      $end  = $arr[ 5 ];
      $bin_ = int( $ini / $bin_size );
      push @{ $genes{$seq_id}{$bin_} }, "$seq_id\t$ini\t$end\tgene";
    }
    close FH;
  }
}

sub loadRegions
{

  # load region coordinates
  warn "loading regions\n" if ( defined $verbose );
  foreach my $file ( @region )
  {
    my $fileh = defineFH( $file );
    open FH, "$fileh" or die "cannot open $file\n";
    while ( <FH> )
    {
      ( $seq_id, $ini, $end ) = split( /\s+/, $_ );
      next unless ( defined $seq{$seq_id} );
      push( @{ $region{$seq_id} }, "$ini-$end" );
    }
    close FH;
  }
}

sub profileSeqs
{

  # compute Kmer frequencies and GC transitions observed
  warn "profiling sequences, kmer=$kmer, window=$win\n" if ( defined $verbose );
  my $bp_slices = 0;
  my @kmer      = createKmer( $kmer - 1, @dna );
  my %kmer      = ();
  my %gct       = ();
  my ( $ini, $end, $reg );
  loadRegions() if ( defined $region );

  while ( ( $seq_id, $seq ) = each %seq )
  {
    warn "    analyzing sequence $seq_id\n" if ( defined $verbose );
    my $last_gc = undef;
    my $len     = length $seq;
    my @slices  = ();

    if ( defined $region )
    {
      next unless ( defined $region{$seq_id} );
      @slices = @{ $region{$seq_id} };
    } else
    {
      push( @slices, "0-$len" );
    }

    foreach $reg ( @slices )
    {
      ( $ini, $end ) = split( /-/, $reg );

      # GC profile
      for ( my $i = $ini ; $i <= $end - $win ; $i += $win )
      {
        $ss = substr( $seq, $i, $win );

        # GC in the window and transition
        my $gc = getBinGC( $seq_id, int( $i + ( $win / 2 ) ) );
        $gct{$last_gc}{$gc}++ if ( defined $last_gc );
        $last_gc = $gc;

        # Effective bases count
        my $nbas = $ss =~ tr/ACGT/ACGT/;
        my $nrep = $ss =~ tr/Rr/Rr/;
        my $nsim = $ss =~ tr/Ss/Ss/;
        my $ntot = $nbas + $nrep + $nsim;
        next unless ( $ntot > 0 );
        push( @{ $ebases{$gc} }, int( 100 * $nbas / $ntot ) );
        push( @{ $rbases{$gc} }, int( 100 * $nrep / $ntot ) );
        push( @{ $sbases{$gc} }, int( 100 * $nsim / $ntot ) );
      }

      # Kmer counts
      for ( my $j = $ini ; $j <= $end - $kmer ; $j++ )
      {
        my $word = substr( $seq, $j, $kmer );
        next if ( $word =~ m/[^ACGT]/ );
        $word = checkRevComp( $word ) if ( defined $revcomp_kmer );
        my $bingc = getBinGC( $seq_id, $j );
        $kmer{$bingc}{$word}++;
      }
    }
  }

  my $kmer_file = "$model.kmer.K$kmer.W$win.data";
  warn "writing k-mer profile in \"$kmer_file\"\n" if ( defined $verbose );
  open K, ">$kmer_file" or die "cannot write file $kmer_file\n";

  # First pass to compute total in GC bins and global total of kmers
  my $tot_sum = 0;
  my %gc_sum  = ();
  foreach my $gc ( @gc )
  {
    foreach my $word ( keys %{ $kmer{$gc} } )
    {
      $tot_sum     += $kmer{$gc}{$word};
      $gc_sum{$gc} += $kmer{$gc}{$word};
    }
  }

  foreach my $gc ( @gc )
  {
    print K "#GC=$gc\n";
    foreach my $w ( @kmer )
    {
      my $tot = 0;
      foreach my $b ( @dna )
      {
        $tot += $kmer{$gc}{"$w$b"} if ( defined $kmer{$gc}{"$w$b"} );
      }
      if ( $tot > 0 )
      {
        foreach my $b ( @dna )
        {
          my $cnt = 0;
          $cnt = $kmer{$gc}{"$w$b"} if ( defined $kmer{$gc}{"$w$b"} );
          my $frq = sprintf( "%.8f", $cnt / $tot );
          my $gc_frq = '-';
          $gc_frq = sprintf( "%.8f", $cnt / $gc_sum{$gc} )
              if ( $gc_sum{$gc} > 0 );

          print K "$w$b    $frq    $gc_frq    $cnt\n";
        }
      } else
      {
        my ( $g, $c ) = split( /-/, $gc );
        my $p = $c / 100;
        my $q = 1 - $p;
        foreach my $b ( @dna )
        {
          my $frq = sprintf( "%.8f", $p / 2 );
          $frq = sprintf( "%.8f", $q / 2 ) if ( $b =~ m/[AT]/ );
          print K "$w$b\t$frq\t0\t0\n";
        }
      }
    }
  }

  my $gc_file = "$model.GCt.W$win.data";
  warn "writing GC transitions in \"$gc_file\"\n" if ( defined $verbose );
  open G, ">$gc_file" or die "cannot open $gc_file\n";
  foreach my $gc1 ( @gc )
  {
    my $tot = 0;
    foreach my $gc2 ( @gc )
    {
      $tot += $gct{$gc1}{$gc2} if ( defined $gct{$gc1}{$gc2} );
    }
    if ( $tot > 0 )
    {
      foreach my $gc2 ( @gc )
      {
        my $cnt = 0;
        $cnt = $gct{$gc1}{$gc2} if ( defined $gct{$gc1}{$gc2} );
        my $frq = sprintf( "%.8f", $cnt / $tot );
        print G "$gc1\t$gc2\t$frq\t$cnt\n";
      }
    } else
    {
      foreach my $gc2 ( @gc )
      {
        my $cnt = 0;
        my $frq = sprintf( "%.8f", 1 / @gc );
        print G "$gc1\t$gc2\t$frq\t$cnt\n";
      }
    }
  }
  close G;
  close K;
}

sub createKmer
{

  # fill a vector with all k-mer combinations
  my $k = shift @_;
  $k--;
  my @old = @_;
  my @new = ();
  if ( $k < 1 )
  {
    return @old;
  } else
  {
    foreach my $e ( @old )
    {
      foreach my $n ( @dna )
      {
        push @new, "$e$n";    # add new element
      }
    }
    createKmer( $k, @new );    # recursive!
  }
}

sub profileRepeats
{

  # call the repeat parsers and write the final table
  warn "profiling repeats\n" if ( defined $verbose );
  loadGenesBin()             if ( defined $gene );
  profileTRF()               if ( defined $trf );
  profileRMAlign()           if ( defined $repeat );

  my $repfile = "$model.repeats.W$win.data";
  warn "writing repeats info in \"$repfile\"\n";
  open R, ">$repfile" or die "cannot open $repfile\n";
  foreach my $gc ( @gc )
  {

    #my $edist = calcPercDist(@{ $ebases{$gc} });
    #my $rdist = calcPercDist(@{ $rbases{$gc} });
    #my $sdist = calcPercDist(@{ $sbases{$gc} });
    my $rep = parseRepeats( @{ $repeat{$gc} } );

  #print R "#GC=$gc    BASES=$edist    REPEATS=$rdist    SIMPLE=$sdist\n$rep\n";
    print R "#GC=$gc\n$rep\n";
  }
  close R;

  my $insfile = "$model.inserts.W$win.data";
  warn "writing repeats inserts info in \"$insfile\"\n";
  open I, ">$insfile" or die "cannot open $insfile\n";
  foreach my $gc ( @gc )
  {
    print I "#GC=$gc\n";
    foreach my $rep1 ( sort keys %{ $repinsert{$gc} } )
    {
      my $sum = $suminsert{$gc}{$rep1};
      foreach my $rep2 ( sort keys %{ $repinsert{$gc}{$rep1} } )
      {
        my $cnt = $repinsert{$gc}{$rep1}{$rep2};
        my $frq = sprintf( "%.8f", $cnt / $sum );
        print I "$rep1\t$rep2\t$frq\t$cnt\n";
      }
    }
  }
  close I;
}

sub parseRepeats
{

  # parse the selected repeats, returns the processed list
  my $res = undef;
  my %rep = ();
  my (
       $rep, $rid, $rfam,  $con, $per, $dir, $div,
       $ins, $del, $indel, $ini, $end, $nfrg
  );
  foreach $rep ( @_ )
  {
    if ( $rep =~ m/SIMPLE/ )
    {
      $res .= "$rep\n";
    } else
    {    # interspersed repeat
      $nfrg = 1;
      if ( $rep =~ m/;/ )
      {
        my @frg = split( /;/, $rep );
        my $frg = shift @frg;
        ( $rid, $rfam, $dir, $div, $ins, $del, $ini, $end ) =
            split( /:/, $frg );
        foreach $frg ( @frg )
        {
          my ( $fdiv, $fins, $fdel, $fini, $fend ) = split( /:/, $frg );
          $div = $fdiv if ( $fdiv > $div );
          $ins = $fins if ( $fins > $ins );
          $del = $fdel if ( $fdel > $del );
          $ini = $fini if ( $fini < $ini );
          $end = $fend if ( $fend > $end );
          $nfrg++;
        }
      } else
      {
        ( $rid, $rfam, $dir, $div, $ins, $del, $ini, $end ) =
            split( /:/, $rep );
      }
      $len = $end - $ini;
      $res .= join ":", $rid, $rfam, $dir, $div, $ins, $del, $len, $nfrg;
      $res .= "\n";
    }
  }
  return $res;
}

sub calcRepDist
{

  # parse the selected repeats, returns the processed list
  my $res  = undef;
  my %rep  = ();
  my %keep = ();
  my $tot  = 0;
  my $fil  = 0;
  my (
       $rep, $rid, $rfam,  $con, $per, $dir, $div,
       $ins, $del, $indel, $ini, $end, $nfrg
  );
  foreach $rep ( @_ )
  {
    $tot++;
    if ( $rep =~ m/SIMPLE/ )
    {
      ( $rfam, $con, $dir, $per, $div, $indel ) = split( /:/, $rep );
      $rep{"$rfam:$con"}{'num'}++;
      $rep{"$rfam:$con"}{'dir'}   .= "$dir,";
      $rep{"$rfam:$con"}{'per'}   .= "$per,";
      $rep{"$rfam:$con"}{'div'}   .= "$div,";
      $rep{"$rfam:$con"}{'indel'} .= "$indel,";
    } else
    {    # interspersed repeat
      $nfrg = 1;
      if ( $rep =~ m/;/ )
      {
        my @frg = split( /;/, $rep );
        my $frg = shift @frg;
        ( $rid, $rfam, $dir, $div, $ins, $del, $ini, $end ) =
            split( /:/, $frg );
        foreach $frg ( @frg )
        {
          my ( $fdiv, $fins, $fdel, $fini, $fend ) = split( /:/, $frg );
          $div = $fdiv if ( $fdiv > $div );
          $ins = $fins if ( $fins > $ins );
          $del = $fdel if ( $fdel > $del );
          $ini = $fini if ( $fini < $ini );
          $end = $fend if ( $fend > $end );
          $nfrg++;
        }
      } else
      {
        ( $rid, $rfam, $dir, $div, $ins, $del, $ini, $end ) =
            split( /:/, $rep );
      }
      my $len = $end - $ini;
      $indel = $ins + $del;
      $rep{"$rid:$rfam"}{'num'}++;
      $rep{"$rid:$rfam"}{'dir'}   .= "$dir,";
      $rep{"$rid:$rfam"}{'div'}   .= "$div,";
      $rep{"$rid:$rfam"}{'indel'} .= "$indel,";
      $rep{"$rid:$rfam"}{'len'}   .= "$len,";
    }
  }

  # create the label with data distributions
  foreach $rep ( keys %rep )
  {
    if ( $rep{$rep}{'num'} < 3 )
    {
      next;
    }
    my @feat = ();
    my %feat = ();

    if ( $rep =~ m/SIMPLE/ )
    {
      @feat = qw/per div indel/;
      my %data = ();
      foreach my $feat ( @feat )
      {
        my @arr = split( /,/, $rep{$rep}{$feat} );
        pop @arr;
        my @arr_sort = sort { $a <=> $b } @arr;
        my ( $q1, $q2, $q3 ) = calcQuartiles( @arr_sort );
        foreach my $x ( @arr )
        {
          my $v = undef;
          if    ( $x <= $q1 ) { $v = "$arr[0]-$q1"; }
          elsif ( $x <= $q2 ) { $v = "$q1-$q2"; }
          elsif ( $x <= $q3 ) { $v = "$q2-$q3"; }
          else { $v = "$q3-$arr[-1]"; }
          push @{ $feat{$feat} }, $v;
        }
      }
      my $n = length scalar @{ $feat{ $feat[ 0 ] } };
      for ( my $i = 0 ; $i <= $n ; $i++ )
      {
        my $per   = $feat{'per'}[ $i ];
        my $div   = $feat{'div'}[ $i ];
        my $indel = $feat{'indel'}[ $i ];
        $data{"$per:$div:$indel"}++;
      }
      foreach my $data ( keys %data )
      {
        if ( $data{$data} > 2 )
        {
          $keep{$rep}{$data} = $data{$data};
          $fil += $data{$data};
        }
      }
    } else
    {
      @feat = qw/div indel len/;
      my %data = ();
      foreach my $feat ( @feat )
      {
        my @arr = split( /,/, $rep{$rep}{$feat} );
        pop @arr;
        my @arr_sort = sort { $a <=> $b } @arr;
        my ( $q1, $q2, $q3 ) = calcQuartiles( @arr_sort );
        foreach my $x ( @arr )
        {
          my $v = undef;
          if    ( $x <= $q1 ) { $v = "$arr[0]-$q1"; }
          elsif ( $x <= $q2 ) { $v = "$q1-$q2"; }
          elsif ( $x <= $q3 ) { $v = "$q2-$q3"; }
          else { $v = "$q3-$arr[-1]"; }
          push @{ $feat{$feat} }, $v;
        }
      }
      my $n = length scalar @{ $feat{ $feat[ 0 ] } };
      for ( my $i = 0 ; $i <= $n ; $i++ )
      {
        my $div = $feat{'div'}[ $i ];
        my $ins = $feat{'indel'}[ $i ];
        my $len = $feat{'len'}[ $i ];
        $data{"$div:$indel:$len"}++;
      }
      foreach my $data ( keys %data )
      {
        if ( $data{$data} > 2 )
        {
          $keep{$rep}{$data} = $data{$data};
          $fil += $data{$data};
        }
      }
    }

    foreach my $rep ( keys %keep )
    {
      foreach my $grp ( keys %{ $keep{$rep} } )
      {
        my $freq = sprintf( "%.10f", $keep{$rep}{$grp} / $fil );
        $res .= "$rep:$freq$grp\n";
      }
    }
  }
  return $res;
}

sub calcBinDist
{

  # define the distribution of values, returns: "min-max"
  my $res = undef;
  my @data = sort { $a <=> $b } @_;
  if ( $data[ 0 ] eq $data[ -1 ] )
  {
    $res = $data[ 0 ];
  } else
  {
    $res = "$data[0]-$data[-1]";
  }
  return $res;
}

sub calcQuartiles
{

  # calculate the quartiles, returns: Q1, Q2, Q3
  my $p2 = int( $#_ / 2 );
  my $p1 = int( $p2 / 2 );
  my $p3 = int( $p2 / 2 ) + $p2;
  return ( $_[ $p1 ], $_[ $p2 ], $_[ $p3 ] );
}

sub profileTRF
{
  # parse TRF output, remove overlapping repeats
  warn "  parsing TRF files\n" if ( defined $verbose );
  my $count = 0;
  foreach my $file ( @trf )
  {
    if ( defined $exclude )
    {
      next if ( $file =~ m/$exclude/ );
    }
    my $last_ini = -1;
    my $last_end = -1;
    my $last_con = 'x';
    my $last_gc  = -1;
    open T, "$file" or die "cannot open $file\n";
    while ( <T> )
    {
      chomp;
      my @line      = split( /\s+/, $_ );
      my $seq_id    = $line[ 0 ];
      my $ini       = $line[ 1 ];
      my $end       = $line[ 2 ];
      if ( $line[3] ne "trf" )
      {
	die "Simple repeat file $file is not in the expected format!\n";
      }
      my $period    = $line[ 5 ];
      my $div       = 100 - $line[ 7 ];
      my $indel     = $line[ 8 ];
      my $consensus = $line[ -1 ];
      my $dir       = '+';
      if ( ! defined $seq{$seq_id} )
      {
	warn "Simple repeat on sequence \'$seq_id\' was not found in genome!\n";
	next;
      }
      my $alt_con = checkRevComp( $consensus );
      if ( $consensus ne $alt_con )
      {
        $consensus = $alt_con;
        $dir       = '-';
      }
      my $label = "SIMPLE:$consensus:$dir:$period:$div:$indel";
      unless ( defined $no_intron )
      {
        next if ( checkGene( $seq_id, $ini, $end ) );
      }

      # Check for overlaping repeats
      if ( $ini >= $last_ini and $ini <= $last_end )
      {
        if ( length $consensus < length $last_con )
        {
          next;    # because last repeat is shorter
        } else
        {
          $count--;
          pop @{ $repeat{$last_gc} };    # last repeat removal
        }
      }

      my $gc = getBinGC( $seq_id, $ini );
      push @{ $repeat{$gc} }, $label;
      $count++;
      $last_ini = $ini;
      $last_end = $end;
      $last_con = $consensus;
      $last_gc  = $gc;
    }
    close T;
  }
  warn "  Profiled $count simple repeats\n" if ( defined $verbose );
}

sub profileRMAlign
{

  # parse RepeatMasker output, mix spliced repeats
  warn "  parsing RepeatMasker align file\n" if ( defined $verbose );
  my %repdata = ();
  my $nf      = 0;
  foreach my $file ( @repeat )
  {
    if ( defined $exclude )
    {
      next if ( $file =~ m/$exclude/ );
    }
    $nf++;
    open T, "$file" or die "cannot open $file\n";
    my @lineRepName = ();
    my %idHash      = ();
    my $lineNum     = 0;
    while ( <T> )
    {

# Must be in this format
# 239 29.42 1.92 0.97 chr1 11678 11780 (249238841) C MER5B#DNA/hAT-Charlie (74) 104 1 m_b1s502i1 4
      $lineNum++;
      chomp;
      next unless ( m/^\d+/ );
      next if ( m/Simple_repeat|Low_complexity|Unknown|Satellite/ );
      my @line   = split( /\s+/, $_ );
      my $div    = $line[ 1 ];
      my $ins    = $line[ 2 ];
      my $del    = $line[ 3 ];
      my $seq_id = $line[ 4 ];
      my $ini    = $line[ 5 ];
      my $end    = $line[ 6 ];

      my ( $dir, $type, $fam, $rini, $rend );
      if ( $line[ 8 ] eq "C" )
      {

        # Negative strand hit
        $dir = "-";
        if ( $line[ 9 ] =~ /(\S+)\#(\S+)/ )
        {
          $type = $1;
          $fam  = $2;
        }
        $rini = $line[ 12 ];
        $rend = $line[ 11 ];
      } else
      {

        # Plus strand hit
        $dir = "+";
        if ( $line[ 8 ] =~ /(\S+)\#(\S+)/ )
        {
          $type = $1;
          $fam  = $2;
        }
        $rini = $line[ 9 ];
        $rend = $line[ 10 ];
      }

      next unless ( defined $seq{$seq_id} );
      my $rid   = "REP.$nf.$line[-1]";
      my $label = "$type:$fam:$dir:$div:$ins:$del:$rini:$rend";

      unless ( defined $no_intron )
      {
        next if ( checkGene( $seq_id, $ini, $end ) );
      }

      ##
      ## Optimsation: track repeat recursions while streaming through
      ## the input file.  Keep a list of ID's associated with line
      ## numbers @lineRepName, and a hash of previous lines where a
      ## repeat ID was found.  Each time we find a pair of joined
      ## fragments we process all the lines in between and add those
      ## to our counts.
      ##
      $lineRepName[ $lineNum ] = "$type#$fam";
      my $id = $line[ $#line ];
      if ( exists $idHash{$id} )
      {
        die "ERROR: $_\n" if ( !exists $repdata{$rid} );
        my $prevLineNum = $idHash{$id};

        # Fix a bug in the previous version which placed the position
        # in the wrong place ( always at the begining of a sequence ).
        my $pos =
            $repdata{$rid}{'ini'} + int( ( $end - $repdata{$rid}{'ini'} ) / 2 );
        my $gc = getBinGC( $repdata{$rid}{'seq_id'}, $pos );
        my $rep_type = $repdata{$rid}{'type'};
        for ( my $i = $prevLineNum + 1 ; $i < $lineNum ; $i++ )
        {
          my $rins_type = $lineRepName[ $i ];
          next if ( !defined $rins_type || $rins_type eq "" );
          $repinsert{$gc}{$rep_type}{$rins_type}++;
          $suminsert{$gc}{$rep_type}++;
        }
      }
      $idHash{$id} = $lineNum;

      if ( defined $repdata{$rid} )
      {
        $repdata{$rid}{'label'} .= ";$div:$ins:$del:$rini:$rend";
        $repdata{$rid}{'end'} = $end;
      } else
      {
        $repdata{$rid}{'label'}  = $label;
        $repdata{$rid}{'ini'}    = $ini;
        $repdata{$rid}{'end'}    = $end;
        $repdata{$rid}{'seq_id'} = $seq_id;
        $repdata{$rid}{'type'}   = "$type#$fam";
      }
    }
    close T;
  }

  foreach my $rid ( keys %repdata )
  {

    # Fix a bug in the previous version which placed the position
    # in the wrong place ( always at the begining of a sequence and more
    # likely to be in the N stretches in an assembly ).
    my $pos = $repdata{$rid}{'ini'} +
        int( ( $repdata{$rid}{'end'} - $repdata{$rid}{'ini'} ) / 2 );
    my $gc = getBinGC( $repdata{$rid}{'seq_id'}, $pos );
    push @{ $repeat{$gc} }, $repdata{$rid}{'label'};
  }
  %repdata = ();
}

sub checkGene
{

  # verify is a region is inside gene annotation
  my ( $seq_id, $ini, $end ) = @_;
  my $res = undef;
  my $seq = substr( $seq{$seq_id}, $ini - 1, $end - $ini );
  $res = 1 if ( $seq =~ m/X/ );
  return $res;
}

sub writeMaskSeq
{

  # write the masked sequence: [ACGT]=effective bases, N=ambiguous bases,
  # S=simple repeats, R=interspersed repeats, X=genes
  my $file = shift @_;
  warn "writing sequence in $file\n" if ( defined $verbose );
  my $good_bases = 0;
  my $null_bases = 0;
  my $rep_bases  = 0;
  my $low_bases  = 0;
  my $func_bases = 0;
  my $tot_bases  = 0;
  open F, ">$file" or die "cannot write $file\n";

  while ( my ( $id, $seq ) = each %seq )
  {
    print F ">$id\n";
    while ( $seq )
    {
      my $s = substr( $seq, 0, 70 );
      print F "$s\n";
      substr( $seq, 0, 70 ) = '';
      $good_bases += $s =~ tr/ACGT//;
      $null_bases += $s =~ tr/N//;
      $rep_bases  += $s =~ tr/R//;
      $low_bases  += $s =~ tr/S//;
      $func_bases += $s =~ tr/X//;
      $tot_bases += length $s;
    }
  }
  close F;
  warn "compressing $file\n" if ( defined $verbose );
  system( "gzip -f $file" );
  warn <<__RES__
    Effective bases            = $good_bases
    Null bases (N)             = $null_bases
    Interspersed repeats bases = $rep_bases
    Low complexity bases       = $low_bases
    Functional bases (genes)   = $func_bases
    Total bases                = $tot_bases
__RES__
      if ( defined $verbose );
}

sub removeTmp
{

  # delete extracted/downloaded files
  warn "removing temporary files\n" if ( defined $verbose );
  my %keep = ();
  if ( defined $keep_dw_files )
  {
    $keep{ $files{$model}{'FAS'} } = 1;
    $keep{ $files{$model}{'RMA'} } = 1;
    $keep{ $files{$model}{'TRF'} } = 1;
    $keep{ $files{$model}{'GEN'} } = 1;
  }

  my @dirs = qw/fasta RM TRF/;

  unlink( $files{$model}{'GEN'} )
      unless ( defined $keep{ $files{$model}{'GEN'} } );
  foreach my $dir ( @dirs )
  {
    chdir $dir;
    opendir D, "." or die "cannot open directory $dir\n";
    while ( my $file = readdir D )
    {
      next if ( defined $keep{$file} );
      next if ( $file eq '.' or $file eq '..' );
      system( "rm -rf $file" );
    }
    closedir D;
    chdir "..";
  }
}

sub calcGC
{

  # calcular GC content (not %)
  my $seq = shift @_;
  my $ngc = $seq =~ tr/CGcg/CGcg/;
  my $nat = $seq =~ tr/ATat/ATat/;
  my $len = $ngc + $nat;
  return 'NA' if ( $len < 1 );
  my $gc = $ngc / $len;
  return $gc;
}

sub classGC
{

  # GC ranges must be equal to @gc values or weird stuff will happen
  my $gc = shift @_;
  return $gc if ( $gc eq 'NA' );
  $gc = sprintf( "%0.3f", $gc );

  my $bin    = 0;
  my $retVal = 'NA';
  foreach my $gcBin ( @gcBins )
  {
    my $gcCut = $gcBin->[ 0 ] / 1000;
    next if ( $gcCut == 0 );
    if ( $gc <= $gcCut )
    {
      $retVal = $gc[ $bin ];
    }
    last if ( $gcCut > $gc );
    $bin++;
  }
  return $retVal;
}

sub calcBinGC
{
  warn "computing GC bins\n" if ( defined $verbose );

  #
  # Categorize the genome into $win (default 1kb) sized
  # bins and count the total of each.
  #
  my @count  = ();
  my %gcHash = ();
  while ( ( $seq_id, $seq ) = each %seq )
  {
    my $len  = length $seq;
    my $gc   = undef;
    my $last = $len - $win;
    for ( my $i = 0 ; $i <= $last ; $i += $win )
    {
      my $s = substr( $seq, $i, $win );
      $gc = calcGC( $s );
      push @{ $gcHash{$seq_id} }, $gc;
      if ( $gc ne "NA" )
      {

        # 0-1000 = 0-100.0%
        $count[ int( 1000 * $gc ) ]++;
      }
    }
    push @{ $gcHash{$seq_id} }, $gc;    # last fragment with length < win
         # RMH: TODO: look into the above statement and fix the last bin
  }

  #
  # Calculate a cumulative 1kb count for each gc level
  #
  my @cumulative;

  # Start with the first count gc=0
  $count[ 0 ] = 0 if ( !defined $count[ 0 ] );
  $cumulative[ 0 ] = $count[ 0 ];
  foreach my $i ( 1 .. $#count )
  {
    $cumulative[ $i ] = 0 if ( !defined $cumulative[ $i ] );
    $count[ $i ]      = 0 if ( !defined $count[ $i ] );
    $cumulative[ $i ] = $cumulative[ $i - 1 ] + $count[ $i ];
  }
  my $total          = $cumulative[ $#cumulative ];
  my $size           = $total / $gcBuckets;
  my $previousbucket = -1;

  #
  #  Create final @gcBins array holding the gc ranges which
  #  evenly divide the genome into $gcBuckets (default 25).
  #
  foreach my $i ( 0 .. $#cumulative )
  {
    my $c      = $cumulative[ $i ];
    my $bucket = int( $c / $size );
    if ( $bucket > $previousbucket )
    {
      push @gcBins, [ $i < $#cumulative ? $i : 1000, $c ];
      $previousbucket = $bucket;
    }
  }

  #
  # setup @gc range strings ( ie. "30.1-35.6", "35.7-41.3" ... )
  #
  @gc = ();
  my $lastCut = -1;
  foreach my $bin ( @gcBins )
  {
    my $gcCut = $bin->[ 0 ] / 10;
    if ( $lastCut >= 0 )
    {
      push @gc, "$lastCut-$gcCut";
    }
    $lastCut = $gcCut;
  }

  #
  # Now calculate the average gc over a larger scale
  # ( binsize = 10kb by default ).
  #
  foreach my $seq_id ( keys( %seq ) )
  {
    my @gc = @{ $gcHash{$seq_id} };
    my $blk = int( ( $binsize / $win ) / 2 );
    for ( my $i = 0 ; $i <= $#gc ; $i++ )
    {
      my $sum = 0;
      my $num = 0;
      my $gc  = 'NA';
      my $ini = $i - $blk;
      $ini = 0 if ( $ini < 0 );
      my $end = $i + $blk;
      $end = $#gc if ( $end > $#gc );
      for ( my $j = $ini ; $j <= $end ; $j++ )
      {
        next unless ( defined $gc[ $j ] );
        next if ( $gc[ $j ] eq 'NA' );
        $sum += $gc[ $j ];
        $num++;
      }
      if ( $num > 0 )
      {
        $gc = classGC( $sum / $num );    # average GC
      }
      push @{ $bingc{$seq_id} }, $gc;
    }
  }

}

sub getBinGC
{

  # returns the GC based on precomputed tables
  my $id  = shift @_;
  my $pos = shift @_;
  my $gc  = $bingc{$id}[ int( $pos / $win ) ];
  return $gc;
}

sub calcPercDist
{

  # returns the distribution of percents
  my $res = undef;
  my @cnt = ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
  my $tot = 0;
  foreach my $x ( @_ )
  {
    $tot++;
    $cnt[ int( $x / 10 ) ]++;
  }
  if ( $tot == 0 )
  {
    return join ",", @cnt;
  }

  # fusion of last range
  my $last = pop @cnt;
  $cnt[ -1 ] += $last;

  foreach my $x ( @cnt )
  {
    my $p = sprintf( "%.4f", $x / $tot );
    $res .= "$p,";
  }
  $res =~ s/,$//;
  return $res;
}

sub calcGCdist
{

  # compute the distribution of GC contents
  my $res  = undef;
  my $tot  = 0;
  my @data = sort { $a <=> $b } @_;
  my ( $q1, $q2, $q3 ) = calcQuartiles( @data );
  my $s1 = 0;
  my $s2 = 0;
  my $s3 = 0;
  my $s4 = 0;
  foreach my $x ( @_ )
  {
    $tot++;
    if    ( $x < $q1 ) { $s1++; }
    elsif ( $x < $q2 ) { $s2++; }
    elsif ( $x < $q3 ) { $s3++; }
    else { $s4++; }
  }
  return 'NA' if ( $tot < 1 );
  $s1 = sprintf( "%.6f", $s1 / $tot );
  $s2 = sprintf( "%.6f", $s2 / $tot );
  $s3 = sprintf( "%.6f", $s3 / $tot );
  $s4 = sprintf( "%.6f", $s4 / $tot );
  $res = "$data[0]-$q1=$s1,$q1-$q2=$s2,$q2-$q3=$s3,$q3-$data[-1]=$s4";
  return $res;
}

sub checkRevComp
{

  # compare a sequence with the reverse/complement,
  # returns the first in alphabetic order
  my $w = shift @_;
  my $r = reverse $w;
  $r =~ tr/ACGTacgt/TGCAtgca/;
  my @w = ( $w, $r );
  @w = sort @w;
  $w = shift @w;
  return $w;
}

sub loadFiles
{

  # UCSC GB files, hash structure is: $files{MODEL}{TYPE} = FILE
  # TYPE can be FAS=Fasta sequences, RMA=RepeatMasker align,
  # TRF=Tandem Repeat Masker out, GEN=Gene annotation
  # Human
  $files{'hg19'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'hg19'}{'RMA'} = 'chromOut.tar.gz';
  $files{'hg19'}{'RMA'} = 'hg19/RepeatMasker-rm405-db20140131/hg19.fa.align.gz';
  $files{'hg19'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'hg19'}{'GEN'} = 'ensGene.txt.gz';

  $files{'hg18'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'hg18'}{'RMA'} = 'chromOut.tar.gz';
  $files{'hg18'}{'RMA'} = 'hg18/RepeatMasker-rm319-db20071204/hg18.fa.align.gz';
  $files{'hg18'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'hg18'}{'GEN'} = 'ensGene.txt.gz';

  $files{'hg17'}{'FAS'} = 'chromFa.zip';
  #$files{'hg17'}{'RMA'} = 'chromOut.zip';
  $files{'hg17'}{'RMA'} = 'hg17/RepeatMasker-rm300-db20040702/hg17.fa.align.gz';
  $files{'hg17'}{'TRF'} = 'chromTrf.zip';
  $files{'hg17'}{'GEN'} = 'ensGene.txt.gz';

  # Cat
  $files{'felCat5'}{'FAS'} = 'felCat5.fa.gz';
  #$files{'felCat5'}{'RMA'} = 'felCat5.fa.out.gz';
  $files{'felCat5'}{'RMA'} = 'felCat5/RepeatMasker-rm405-db20140131/felCat5.fa.align.gz';
  $files{'felCat5'}{'TRF'} = 'felCat5.trf.bed.gz';
  $files{'felCat5'}{'GEN'} = 'refGene.txt.gz';

  # Chicken
  $files{'galGal3'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'galGal3'}{'RMA'} = 'chromOut.tar.gz';
  $files{"galGal3"}{"RMA"} = "galGal3/RepeatMasker-rm327-db20090204/galGal3.fa.align.gz";
  $files{'galGal3'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'galGal3'}{'GEN'} = 'ensGene.txt.gz';

  # Chimpanzee
  $files{'panTro4'}{'FAS'} = 'panTro4.fa.gz';
  #$files{'panTro4'}{'RMA'} = 'panTro4.fa.out.gz';
  $files{"panTro4"}{"RMA"} = "panTro4/RepeatMasker-rm405-db20140131/panTro4.fa.align.gz";
  $files{'panTro4'}{'TRF'} = 'panTro4.trf.bed.gz';
  $files{'panTro4'}{'GEN'} = 'refGene.txt.gz';

  # Cow
  $files{'bosTau4'}{'FAS'} = 'bosTau4.fa.gz';
  #$files{'bosTau4'}{'RMA'} = 'bosTau4.fa.out.gz';
  $files{"bosTau4"}{"RMA"} = "bosTau4/RepeatMasker-rm328-db20090604/bosTau4.fa.align.gz";
  $files{'bosTau4'}{'TRF'} = 'bosTau4.trf.bed.gz';
  $files{'bosTau4'}{'GEN'} = 'ensGene.txt.gz';

  # Dog
  $files{'canFam2'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'canFam2'}{'RMA'} = 'chromOut.tar.gz';
  $files{"canFam2"}{"RMA"} = "canFam2/RepeatMasker-rm330-db20120124/canFam2.fa.align.gz";
  $files{'canFam2'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'canFam2'}{'GEN'} = 'ensGene.txt.gz';

  # Elephant
  $files{'loxAfr3'}{'FAS'} = 'loxAfr3.fa.gz';
  #$files{'loxAfr3'}{'RMA'} = 'loxAfr3.fa.out.gz';
  $files{"loxAfr3"}{"RMA"} = "loxAfr3/RepeatMasker-rm405-db20140131/loxAfr3.fa.align.gz";
  $files{'loxAfr3'}{'TRF'} = 'loxAfr3.trf.bed.gz';
  $files{'loxAfr3'}{'GEN'} = 'ensGene.txt.gz';

  # Fugu
  $files{'fr2'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'fr2'}{'RMA'} = 'chromOut.tar.gz';
  $files{"fr2"}{"RMA"} = "fr2/RepeatMasker-rm325-db20080611/fr2.fa.align.gz";
  $files{'fr2'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'fr2'}{'GEN'} = 'ensGene.txt.gz';

  # Guinea Pig
  $files{'cavPor3'}{'FAS'} = 'cavPor3.fa.gz';
  #$files{'cavPor3'}{'RMA'} = 'cavPor3.fa.out.gz';
  $files{"cavPor3"}{"RMA"} = "cavPor3/RepeatMasker-rm405-db20140131/cavPor3.fa.align.gz";
  $files{'cavPor3'}{'TRF'} = 'cavPor3.trf.bed.gz';
  $files{'cavPor3'}{'GEN'} = 'ensGene.txt.gz';

  # Horse
  $files{'equCab2'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'equCab2'}{'RMA'} = 'chromOut.tar.gz';
  $files{"equCab2"}{"RMA"} = "equCab2/RepeatMasker-rm405-db20140131/equCab2.fa.align.gz";
  $files{'equCab2'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'equCab2'}{'GEN'} = 'ensGene.txt.gz';

  # Lizard
  $files{'anoCar2'}{'FAS'} = 'anoCar2.fa.gz';
  #$files{'anoCar2'}{'RMA'} = 'anoCar2.fa.out.gz';
  $files{"anoCar2"}{"RMA"} = "anoCar2/RepeatMasker-rm405-db20140131/anoCar2.fa.align.gz";
  $files{'anoCar2'}{'TRF'} = 'anoCar2.trf.bed.gz';
  $files{'anoCar2'}{'GEN'} = 'ensGene.txt.gz';

  # Marmoset
  $files{'calJac3'}{'FAS'} = 'calJac3.fa.gz';
  #$files{'calJac3'}{'RMA'} = 'calJac3.fa.out.gz';
  $files{"calJac3"}{"RMA"} = "calJac3/RepeatMasker-rm405-db20140131/calJac3.fa.align.gz";
  $files{'calJac3'}{'TRF'} = 'calJac3.trf.bed.gz';
  $files{'calJac3'}{'GEN'} = 'ensGene.txt.gz';

  # Medaka
  #$files{'oryLat2'}{'FAS'} = 'oryLat2.fa.gz';
  #$files{'oryLat2'}{'RMA'} = 'oryLat2.fa.out.gz';
  #$files{'oryLat2'}{'TRF'} = 'oryLat2.trf.bed.gz';
  #$files{'oryLat2'}{'GEN'} = 'ensGene.txt.gz';

  # Mouse
  $files{'mm10'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'mm10'}{'RMA'} = 'chromOut.tar.gz';
  $files{"mm10"}{"RMA"} = "mm10/RepeatMasker-rm405-db20140131/mm10.fa.align.gz";
  $files{'mm10'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'mm10'}{'GEN'} = 'ensGene.txt.gz';

  $files{'mm9'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'mm9'}{'RMA'} = 'chromOut.tar.gz';
  $files{"mm9"}{"RMA"} = "mm9/RepeatMasker-rm328-db20090604/mm9.fa.align.gz";
  $files{'mm9'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'mm9'}{'GEN'} = 'ensGene.txt.gz';

  # Oposum
  $files{'monDom5'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'monDom5'}{'RMA'} = 'chromOut.tar.gz';
  $files{"monDom5"}{"RMA"} = "monDom5/RepeatMasker-rm405-db20140131/monDom5.fa.align.gz";
  $files{'monDom5'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'monDom5'}{'GEN'} = 'ensGene.txt.gz';

  # Orangutan
  $files{'ponAbe2'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'ponAbe2'}{'RMA'} = 'chromOut.tar.gz';
  $files{"ponAbe2"}{"RMA"} = "ponAbe2/RepeatMasker-rm405-db20140131/ponAbe2.fa.align.gz";
  $files{'ponAbe2'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'ponAbe2'}{'GEN'} = 'ensGene.txt.gz';

  # Panda
  $files{'ailMel1'}{'FAS'} = 'ailMel1.fa.gz';
  #$files{'ailMel1'}{'RMA'} = 'ailMel1.fa.out.gz';
  $files{"ailMel1"}{"RMA"} = "ailMel1/RepeatMasker-rm405-db20140131/ailMel1.fa.align.gz";
  $files{'ailMel1'}{'TRF'} = 'ailMel1.trf.bed.gz';
  $files{'ailMel1'}{'GEN'} = 'ensGene.txt.gz';

  # Pig
  $files{'susScr2'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'susScr2'}{'RMA'} = 'chromOut.tar.gz';
  $files{"susScr2"}{"RMA"} = "susScr2/RepeatMasker-rm330-db20120124/susScr2.fa.align.gz";
  $files{'susScr2'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'susScr2'}{'GEN'} = 'ensGene.txt.gz';

  # Platypus
  $files{'ornAna1'}{'FAS'} = 'ornAna1.fa.gz';
  #$files{'ornAna1'}{'RMA'} = 'ornAna1.fa.out.gz';
  $files{"ornAna1"}{"RMA"} = "ornAna1/RepeatMasker-rm405-db20140131/ornAna1.fa.align.gz";
  $files{'ornAna1'}{'TRF'} = 'ornAna1.trf.bed.gz';
  $files{'ornAna1'}{'GEN'} = 'ensGene.txt.gz';

  # Rabbit
  $files{'oryCun2'}{'FAS'} = 'oryCun2.fa.gz';
  #$files{'oryCun2'}{'RMA'} = 'oryCun2.fa.out.gz';
  $files{"oryCun2"}{"RMA"} = "oryCun2/RepeatMasker-rm405-db20140131/oryCun2.fa.align.gz";
  $files{'oryCun2'}{'TRF'} = 'oryCun2.trf.bed.gz';
  $files{'oryCun2'}{'GEN'} = 'ensGene.txt.gz';

  # Rat
  $files{'rn4'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'rn4'}{'RMA'} = 'chromOut.tar.gz';
  $files{"rn4"}{"RMA"} = "rn4/RepeatMasker-rm330-db20120124/rn4.fa.align.gz";
  $files{'rn4'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'rn4'}{'GEN'} = 'ensGene.txt.gz';

  # Rhesus
  $files{'rheMac2'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'rheMac2'}{'RMA'} = 'chromOut.tar.gz';
  $files{"rheMac2"}{"RMA"} = "rheMac2/RepeatMasker-rm319-db20071204/rheMac2.fa.align.gz";
  $files{'rheMac2'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'rheMac2'}{'GEN'} = 'ensGene.txt.gz';

  # Sheep
  #$files{'oviAri1'}{'FAS'} = 'oviAri1.fa.gz';
  #$files{'oviAri1'}{'RMA'} = 'oviAri1.fa.out.gz';
  #$files{'oviAri1'}{'TRF'} = 'oviAri1.trf.bed.gz';
  #$files{'oviAri1'}{'GEN'} = 'refGene.txt.gz';

  # Stickleback
  $files{'gasAcu1'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'gasAcu1'}{'RMA'} = 'chromOut.tar.gz';
  $files{"gasAcu1"}{"RMA"} = "gasAcu1/RepeatMasker-rm405-db20140131/gasAcu1.fa.align.gz";
  $files{'gasAcu1'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'gasAcu1'}{'GEN'} = 'ensGene.txt.gz';

  # Tetraodon
  #$files{'tetNig2'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'tetNig2'}{'RMA'} = 'chromOut.tar.gz';
  #$files{'tetNig2'}{'TRF'} = 'chromTrf.tar.gz';
  #$files{'tetNig2'}{'GEN'} = 'ensGene.txt.gz';

  # X. tropicalis
  $files{'xenTro2'}{'FAS'} = 'xenTro2.fa.gz';
  #$files{'xenTro2'}{'RMA'} = 'xenTro2.rmsk.out.gz';
  $files{"xenTro2"}{"RMA"} = "xenTro2/RepeatMasker-rm327-db20090202/xenTro2.fa.align.gz";
  $files{'xenTro2'}{'TRF'} = 'xenTro2.trf.bed.gz';
  $files{'xenTro2'}{'GEN'} = 'ensGene.txt.gz';

  # Zebra finch
  $files{'taeGut1'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'taeGut1'}{'RMA'} = 'chromOut.tar.gz';
  $files{"taeGut1"}{"RMA"} = "taeGut1/RepeatMasker-rm405-db20140131/taeGut1.fa.align.gz";
  $files{'taeGut1'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'taeGut1'}{'GEN'} = 'ensGene.txt.gz';

  # Zefrafish
  $files{'danRer7'}{'FAS'} = 'danRer7.fa.gz';
  #$files{'danRer7'}{'RMA'} = 'danRer7.fa.out.gz';
  $files{"danRer7"}{"RMA"} = "danRer7/RepeatMasker-rm405-db20140131/danRer7.fa.align.gz";
  $files{'danRer7'}{'TRF'} = 'danRer7.trf.bed.gz';
  $files{'danRer7'}{'GEN'} = 'ensGene.txt.gz';

  # C. intestinalis
  $files{'ci2'}{'FAS'} = 'ScaffoldFa.zip';
  #$files{'ci2'}{'RMA'} = 'Scaffold.out.zip';
  $files{"ci2"}{"RMA"} = "ci2/RepeatMasker-rm405-db20140131/ci2.fa.align.gz";
  $files{'ci2'}{'TRF'} = 'ScaffoldTrf.zip';
  $files{'ci2'}{'GEN'} = 'ensGene.txt.gz';

  # S. purpuratus
  $files{'strPur2'}{'FAS'} = 'strPur2.fa.gz';
  #$files{'strPur2'}{'RMA'} = 'strPur2.fa.out.gz';
  $files{"strPur2"}{"RMA"} = "strPur2/RepeatMasker-rm405-db20140131/strPur2.fa.align.gz";
  $files{'strPur2'}{'TRF'} = 'strPur2.trf.bed.gz';
  $files{'strPur2'}{'GEN'} = 'refGene.txt.gz';

  # A. gambiae
  $files{'anoGam1'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'anoGam1'}{'RMA'} = 'chromOut.tar.gz';
  $files{"anoGam1"}{"RMA"} = "anoGam1/RepeatMasker-rm405-db20140131/anoGam1.fa.align.gz";
  $files{'anoGam1'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'anoGam1'}{'GEN'} = 'ensGene.txt.gz';

  # A. mellifera
  $files{'apiMel3'}{'FAS'} = 'GroupFa.zip';
  #$files{'apiMel3'}{'RMA'} = 'GroupOut.zip';
  $files{"apiMel3"}{"RMA"} = "apiMel3/RepeatMasker-rm405-db20140131/apiMel3.fa.align.gz";
  $files{'apiMel3'}{'TRF'} = 'GroupTrf.zip';
  $files{'apiMel3'}{'GEN'} = 'ensGene.txt.gz';

  # D. melanogaster
  $files{'dm3'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'dm3'}{'RMA'} = 'chromOut.tar.gz';
  $files{"dm3"}{"RMA"} = "dm3/RepeatMasker-rm405-db20140131/dm3.fa.align.gz";
  $files{'dm3'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'dm3'}{'GEN'} = 'ensGene.txt.gz';

  # C. elegans
  $files{'ce10'}{'FAS'} = 'chromFa.tar.gz';
  #$files{'ce10'}{'RMA'} = 'chromOut.tar.gz';
  $files{"ce10"}{"RMA"} = "ce10/RepeatMasker-rm405-db20140131/ce10.fa.align.gz";
  $files{'ce10'}{'TRF'} = 'chromTrf.tar.gz';
  $files{'ce10'}{'GEN'} = 'ensGene.txt.gz';
}

# EOF
