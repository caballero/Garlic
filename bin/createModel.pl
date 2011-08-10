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
    -w --win           Profile window size           Integer        1000
	-b --binsize       Bin size (for %GC)            Integer        20000
    
    -f --fasta         Fasta sequences               FileName*
    -r --repeats       RepeatMasker output           FileName*
    -t --trf           TRF output                    FileName*
    -g --genes         Gene annotation               FileName*
    -e --exclude       Exclude this sequences        Pattern**
    
    --no_mask_gene     Don't mask genes
    --no_mask_repeat   Don't mask repeats
    --no_mask_trf      Don't mask simple repeats
    --no_repeat_table  Don't create repeats table
    --no_kmer_table    Don't compute k-mer composition
    --rm_tmp           Remove temp files
    --keep_dw_files    Keep downloaded files (overrules --rm_tmp)
    --gc_post_mask     Compute GC bins after masking (default is before)
    --revcomp_kmer     Count kmers in reverse complement chain

  * File can be compressed (.gz/.bz2), if you are passing more than one file 
    separate them with ',' Example: "-f chr1.fa,chr2.fa,chr3.fa"
 ** This pattern will match the fasta files. Example: "-e hap1" will exclude 
    all files with "hap1" in the name.
    
=head1 EXAMPLES

You can obtain the models from our website, we are currently suporting some 
organisms with complete annotation in the UCSC Genome Database:
<http://hgdownload.cse.ucsc.edu/downloads.html>

You can create your own model fetching the data from the UCSC site:

  perl createModel.pl -m hg19

The model names are according to UCSC:
  - hg19    -> human genome release 19 (GRCh37)
  - hg18    -> human genome release 18 (GRCh36)
  - mm9     -> mouse genome release 9
  - equCab2 -> horse genome release 2
  - ...

This script will download the required files and process it. Caution, the files
are really big and the download time may be large.

Alternatively, you can use your own sequences and annotations to create a model:

  perl createModel.pl -m myOrg -f myOrg.fa -r RM.out -t TRF.out -g Genes.table

In this case you must provide the sequences in a Fasta file, the RepeatMasker 
output, the Tandem Repeat Finder output and the annotated genes in a tabular 
text file similar to the knownGenes.txt format.
  
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
my $model           =  undef;      # Model definition
my $dir             = 'data';      # Output directory
my $kmer            =      4;      # K-mer size
my $win             =   1000;      # Window size (non-overlapping)
my $binsize         =   2000;      # Bin size
my $fasta           =  undef;      # Sequences in fasta files
my $repeat          =  undef;      # RepeatMasker output
my $trf             =  undef;      # TRF output
my $gene            =  undef;      # Gene annotation
my $help            =  undef;      # Help flag
my $verbose         =  undef;      # Verbose mode flag
my $rm_tmp          =  undef;      # Remove downloaded files
my $exclude         =  undef;      # Exclude sequence name pattern
my $mask_repeat     =  undef;      # Repeat masking flag
my $mask_trf        =  undef;      # TRF masking flag
my $no_mask_gene    =  undef;      # Gene masking flag
my $no_mask_repeat  =  undef;      # Repeats masking flag
my $no_mask_trf     =  undef;      # TRF masking flag
my $no_repeat_table =  undef;      # No repeats profile flag
my $no_kmer_table   =  undef;      # No kmer profile flag
my $write_mask_seq  =  undef;      # Write fasta flag
my $keep_dw_files   =  undef;      # Keep downloaded files
my $gc_post_mask    =  undef;      # Compute GC bins before masking
my $revcomp_kmer    =  undef;      # Count kmers in the reverse-complement chain

# Fetch options
GetOptions(
    'h|help'             => \$help,
    'v|verbose'          => \$verbose,
    'm|model=s'          => \$model,
    'd|dir:s'            => \$dir,
    'k|kmer:i'           => \$kmer,
    'w|win:i'            => \$win,
	'b|binsize:i'        => \$binsize,
    'f|fasta:s'          => \$fasta,
    'r|repeat:s'         => \$repeat,
    't|trf:s'            => \$trf,
    'g|gene:s'           => \$gene,
    'e|exclude:s'        => \$exclude,
    'rm_tmp'             => \$rm_tmp,
    'no_mask_gene'       => \$no_mask_gene,
    'no_mask_repeat'     => \$no_mask_repeat,
    'no_mask_trf'        => \$no_mask_trf,
    'no_repeat_table'    => \$no_repeat_table,
    'no_kmer_table'      => \$no_kmer_table,
    'write_mask_seq'     => \$write_mask_seq,
    'keep_dw_files'      => \$keep_dw_files,
    'gc_post_mask'       => \$gc_post_mask,
    'revcomp_kmer'       => \$revcomp_kmer
) or pod2usage(-verbose => 2);

# Call help if required
pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $model);

# UCSC model files
my %files = ();
loadFiles();

# Configurable parameters
my $get           = 'wget -c';           # command to fetch files from internet
   $get          .= ' -q' unless (defined $verbose);
my $unpack        = 'tar zxf';           # command to unpack the files downloaded
my $unzip         = 'unzip';             # command to decompress the files downloaded
my $gunzip        = 'gunzip -c';         # command to decompress the files downloaded
my $ucsc          = 'http://hgdownload.cse.ucsc.edu/goldenPath'; # UCSC url
my $ucsc_genome   = "$ucsc/$model/bigZips/"  . $files{$model}{'FAS'};
my $ucsc_repeat   = "$ucsc/$model/bigZips/"  . $files{$model}{'RMO'};
my $ucsc_trf      = "$ucsc/$model/bigZips/"  . $files{$model}{'TRF'};
my $ucsc_gene     = "$ucsc/$model/database/" . $files{$model}{'GEN'}; 

# Main variables
my @fasta    = ();
my @repeat   = ();
my @trf      = ();
my @gene     = ();
my %bingc    = ();
my %seq      = ();
my %repeat   = ();
my %genes    = ();
my @dna      = qw/A C G T/;
my ($seq, $ss, $seq_id, $ini, $end, $len);
# GC content bins, if you change this array, change the classGC() subroutine too.
my @gc = qw/0-5 5-10 10-15 15-20 20-25 25-30 30-35 35-40 40-45 45-50 50-55 55-60 60-65 65-70 70-75 75-80 80-85 85-90 90-95 95-100/;

# Check directories, create them if required
unless (-e "$dir" and -d "$dir") {
    warn "creating directory $dir\n" if (defined $verbose);
    mkdir "$dir";
}
unless (-e "$dir/$model" and -d "$dir/$model") {
    warn "creating directory $dir/$model\n" if (defined $verbose);
    mkdir "$dir/$model";
}
warn "moving to $dir/$model\n" if (defined $verbose);
chdir "$dir/$model" or die "cannot move to $dir/$model";

# Get files from UCSC genome database if required
getUCSC_fasta()  unless (defined $fasta);
getUCSC_repeat() unless (defined $repeat or (defined $no_repeat_table and defined $no_mask_repeat));
getUCSC_trf()    unless (defined $trf    or (defined $no_repeat_table and defined $no_mask_trf));
getUCSC_gene()   unless (defined $gene   or  defined $no_mask_gene);

# Creating list of files to process
@fasta   = split (/,/, $fasta);
@repeat  = split (/,/, $repeat);
@trf     = split (/,/, $trf);
@gene    = split (/,/, $gene);
$exclude =~ s/,/|/g if (defined $exclude);

# Load sequences and mask it
readFasta();
calcBinGC()  unless (defined $gc_post_mask);
maskRepeat() unless (defined $no_mask_repeat);
maskTRF()    unless (defined $no_mask_trf);
maskGene()   unless (defined $no_mask_gene);
calcBinGC()  if     (defined $gc_post_mask);
writeMaskSeq("$model.masked.fa") if (defined $write_mask_seq);

# Create the K-mer/Window and GCt tables
profileSeqs() unless (defined $no_kmer_table);

# Create the Repeat Table
profileRepeats() unless (defined $no_repeat_table);

# Create the model info file
writeModelInfo();

# Clean up
removeTmp() if (defined $rm_tmp);

warn "Done\n" if (defined $verbose);

#################################################
##      S  U  B  R  O  U  T  I  N  E  S        ##
#################################################

sub searchFiles {
    # deep search of files/patterns in a directory
    my ($pat, $dir) = @_;
    my $files = undef;
    find ( sub { $files .= $File::Find::name . ',' if (m/$pat/) }, $dir );
    return $files;
}

sub getUCSC_gene {
    # Gene models to mask functional regions
    $gene = $files{$model}{'GEN'};
    die "Sorry, $model don't have annotation registrated\n" unless (defined $gene);

    if (-e $gene) {
        warn "$gene present, using it\n" if (defined $verbose);
    }
    else {
        warn "obtaining Gene files from $ucsc\n" if (defined $verbose);
        system ("$get $ucsc_gene");
        die "cannot find Gene output in $ucsc_gene" unless (-e $gene and -s $gene);
    }
}

sub getUCSC_trf {
    # TRF output is used for simple repeat annotation
    my $target_file = $files{$model}{'TRF'};
    die "Sorry, $model don't have TRF registrated\n" unless (defined $target_file);

    if (-e "TRF/$target_file") {
        warn "$target_file present, using it\n" if (defined $verbose);
        chdir 'TRF' or die "cannot move to TRF directory\n";

    }
    else {
        warn "obtaining TRF files from $ucsc\n" if (defined $verbose);
        mkdir 'TRF' unless (-e 'TRF' and -d 'TRF');
        chdir 'TRF' or die "cannot move to TRF directory\n";
        system ("$get $ucsc_trf");
        die "cannot find TRF output in $ucsc_trf" unless (-e $target_file);
    }
    
    if ($target_file =~ m/tar.gz$/) {
        warn "unpacking TAR\n" if (defined $verbose);
        system ("$unpack $target_file");
    }
    elsif ($target_file =~ m/.zip$/) {
        warn "unpacking ZIP\n" if (defined $verbose);
        system ("$unzip $target_file");
    }
    elsif($target_file =~ m/.gz$/) {
        warn "unzipping GZ\n" if (defined $verbose);
		my $out_file = $target_file;
		$out_file =~ s/.gz$//;
        system ("$gunzip $target_file > $out_file");
    }
    else {
        die "cannot recognize the file format of $target_file\n";
    }
    
    chdir '..';
    warn "searching TRF files\n" if (defined $verbose);
    $trf = searchFiles('.bed$', 'TRF');
    die "cannot find TRF files!" unless (defined $trf);
}

sub getUCSC_repeat {
    # RepeatMasker output
    my $target_file = $files{$model}{'RMO'};
    die "Sorry, $model don't have RepeatMasker registrated\n" unless (defined $target_file);

    if (-e "RM/$target_file") {
        warn "$target_file present, using it\n" if (defined $verbose);
        chdir 'RM' or die "cannot move to RM directory\n";

    }
    else {
        warn "obtaining RepeatMasker files from $ucsc\n" if (defined $verbose);
        mkdir 'RM' unless (-e 'RM' and -d 'RM');
        chdir 'RM' or die "cannot move to RM directory\n";
        system ("$get $ucsc_repeat");
        die "cannot find RepeatMasker output in $ucsc_repeat" unless (-e $target_file);
    }
    
    if ($target_file =~ m/tar.gz$/) {
        warn "unpacking TAR\n" if (defined $verbose);
        system ("$unpack $target_file");
    }
    elsif ($target_file =~ m/.zip$/) {
        warn "unpacking ZIP\n" if (defined $verbose);
        system ("$unzip $target_file");
    }
    elsif($target_file =~ m/.gz$/) {
        warn "unzipping GZ\n" if (defined $verbose);
		my $out_file = $target_file;
		$out_file =~ s/.gz$//;
        system ("$gunzip $target_file > $out_file");
    }
    else {
        die "cannot recognize the file format of $target_file\n";
    }
    
    chdir '..';
    warn "searching RM files\n" if (defined $verbose);
    $repeat = searchFiles('.out$', 'RM');
    die "cannot find RM files!" unless (defined $repeat);
}

sub getUCSC_fasta {
    # Chromosomal sequences
    my $target_file = $files{$model}{'FAS'};
    die "Sorry, $model don't have Fasta registrated\n" unless (defined $target_file);

    if (-e "fasta/$target_file") {
        warn "$target_file present, using it\n" if (defined $verbose);
        chdir 'fasta' or die "cannot move to fasta directory\n";
    }
    else {
        warn "obtaining fasta files from $ucsc\n" if (defined $verbose);
        mkdir 'fasta' unless (-e 'fasta' and -d 'fasta');
        chdir 'fasta' or die "cannot move to fasta directory\n";
        system ("$get $ucsc_genome");
        die "cannot find genomic sequences in $ucsc_genome" unless (-e $target_file);
    }
    
    if ($target_file =~ m/tar.gz$/) {
        warn "unpacking TAR\n" if (defined $verbose);
        system ("$unpack $target_file");
    }
    elsif ($target_file =~ m/.zip$/) {
        warn "unpacking ZIP\n" if (defined $verbose);
        system ("$unzip $target_file");
    }
    elsif($target_file =~ m/.gz$/) {
        warn "unzipping GZ\n" if (defined $verbose);
		my $out_file = $target_file;
		$out_file =~ s/.gz$//;
        system ("$gunzip $target_file > $out_file");
    }
    else {
        die "cannot recognize the file format of $target_file\n";
    }
    
    chdir '..';
    warn "searching fasta files\n" if (defined $verbose);
    $fasta = searchFiles('.fa$', 'fasta');
    die "cannot find fasta files!" unless (defined $fasta);
}

sub defineFH {
	my $fi = shift @_;
	my $fh = $fi;
    $fh = "gunzip  -c $fi | " if ($fi =~ m/\.gz$/);
    $fh = "bunzip2 -c $fi | " if ($fi =~ m/\.bz2$/);
	return $fh;
}
		
sub readFasta {
    # read tha fasta files, apply filters if required
    warn "loading fasta sequences\n" if (defined $verbose);
    foreach my $file (@fasta) {
        if (defined $exclude) {
            next if ($file =~ m/$exclude/);
        }
        warn "  reading $file\n" if (defined $verbose);
        my $fileh = defineFH($file);
        open FH, "$fileh" or die "cannot open $file\n";
        while (<FH>) {
            chomp;
            if (m/^>(.+)/) { 
                $seq_id = $1;        
            }
            else {
                s/[^ACGTacgt]/N/g;
                $seq{$seq_id} .= $_; 
            }
        }
        close FH;
    }
}

sub maskRepeat {
    # mask sequences with RepeatMasker annotation
    warn "masking repeats\n" if (defined $verbose);
    foreach my $file (@repeat) {
        if (defined $exclude) {
            next if ($file =~ m/$exclude/);
        }
        my $fileh = defineFH($file);
        open FH, "$fileh" or die "cannot open $file\n";
        while (<FH>) {
            s/^\s*//;
            next unless (m/^\d+/); # skip headers
            my @arr = split (/\s+/, $_);
            $seq_id = $arr[4];
            next unless (defined $seq{$seq_id});
            $ini = $arr[5];
            $end = $arr[6];
            $len = $end - $ini + 1;
            substr ($seq{$seq_id}, $ini - 1, $len) = 'R' x $len;
        }
        close FH;
    }
}

sub maskTRF {
    # mask sequences with TRF annotation
    warn "masking simple repeats\n" if (defined $verbose);
    foreach my $file (@trf) {
        if (defined $exclude) {
            next if ($file =~ m/$exclude/);
        }
        my $fileh = defineFH($file);
        open FH, "$fileh" or die "cannot open $file\n";
        while (<FH>) {
            my @arr = split (/\t/, $_);
            $seq_id = $arr[0];
            next unless (defined $seq{$seq_id});
            $ini = $arr[1];
            $end = $arr[2];
            $len = $end - $ini;
            substr ($seq{$seq_id}, $ini, $len) = 'S' x $len;
        }
        close FH;
    }
}

sub maskGene {
    # mask sequences with gene annotation
    warn "masking genes\n" if (defined $verbose);
    foreach my $file (@gene) {
        my $fileh = defineFH($file);
        open FH, "$fileh" or die "cannot open $file\n";
        while (<FH>) {
            my @arr = split (/\t/, $_);
            $seq_id = $arr[2];
            next unless (defined $seq{$seq_id});
            $ini = $arr[4];
            $end = $arr[5];
            $len = $end - $ini - 1;
            substr ($seq{$seq_id}, $ini, $len) = 'X' x $len;
            push @{ $genes{$seq_id} }, "$ini-$end";
        }
        close FH;
    }
}

sub profileSeqs {
    # compute Kmer frequencies and GC transitions observed
    warn "profiling sequences, kmer=$kmer and window=$win\n" if (defined $verbose);
    my $bp_slices = 0;
    my @kmer      = createKmer($kmer - 1, @dna);
    my %kmer      = ();
    my %gct       = ();
    
    foreach $seq_id (keys %seq) {
	    warn "  analyzing sequence $seq_id\n" if (defined $verbose);
        $seq = $seq{$seq_id};
	    $seq =~ s/[^ACGT]/N/g;
	    my $last_gc = undef;
	    my $len = length $seq;
		# GC windows transitions
		my $last = $len - $win;
	    for (my $i = 0; $i <= $last; $i += $win) {
	        $ss = substr ($seq, $i, $win);
	        my $num_N = $ss =~ tr/N/N/;
	        my $frq_N = $num_N / (length $ss);
	        next if ($frq_N > 0.3); # No more than 30% bad bases 
	        
	        # GC in the window and transition
	        my $gc = getBinGC($seq_id, int($i + ($win / 2)) );
	        $gct{$last_gc}{$gc}++ if (defined $last_gc);
	        $last_gc = $gc;   
	    }
		# Kmer counts
		$last = $len - $kmer;
		for (my $j = 0; $j <= $last; $j++) {
	        my $word = substr ($seq, $j, $kmer);
			next if ($word =~ m/[^ACGT]/);
			$word = checkRevComp($word) if (defined $revcomp_kmer);
			my $bingc = getBinGC($seq_id, $j);
			$kmer{$bingc}{$word}++;
	    }
    }
    
    my $kmer_file = "$model.kmer.K$kmer.W$win.data";
    warn "writing k-mer profile in \"$kmer_file\"\n" if (defined $verbose);
    open K, ">$kmer_file" or die "cannot write file $kmer_file\n";
    
    # First pass to compute total in GC bins and global total of kmers
    my $tot_sum = 0;
    my %gc_sum  = ();
    foreach my $gc (@gc) {
        foreach my $word (keys %{ $kmer{$gc} }) {
            $tot_sum     += $kmer{$gc}{$word};
            $gc_sum{$gc} += $kmer{$gc}{$word};
        }
    }
    
    foreach my $gc (@gc) {
	    print K "#GC=$gc\n";
	    foreach my $word (@kmer) {
	        my $tot = 0;
	        foreach my $b (@dna) {
	            $tot += $kmer{$gc}{"$word$b"} if (defined $kmer{$gc}{"$word$b"});
	        }
	        if ($tot > 0) {
	            foreach my $b (@dna) {
	                my $cnt     = 0;
	                   $cnt     = $kmer{$gc}{"$word$b"} if (defined $kmer{$gc}{"$word$b"});
	                my $frq     = sprintf ("%.8f", $cnt / $tot);
	                my $gc_frq  = '-';
	                   $gc_frq  = sprintf ("%.8f", $cnt / $gc_sum{$gc}) if ($gc_sum{$gc} > 0);
	                
	                my $tot_frq = '-';
	                   $tot_frq = sprintf ("%.8f", $cnt / $tot_sum);
	                   
	                print K "$word$b\t$frq\t$gc_frq\t$tot_frq\t$cnt\n";
	            }
	        }
	        else {
				my ($g, $c) = split (/-/, $gc);
				my $p = ($c - 5) / 100;
				my $q = 1 - $p;
	            foreach my $b (@dna) {
	                my $frq = sprintf ("%.8f", $p / 2);
					   $frq = sprintf ("%.8f", $q / 2) if ($b =~ m/[AT]/);
	                print K "$word$b\t$frq\t0\t0\t0\n";   
	            }
	        }
	    }
	}
	
	my $gc_file = "$model.GCt.W$win.data";
	warn "writing GC transitions in \"$gc_file\"\n" if (defined $verbose);
	open G, ">$gc_file" or die "cannot open $gc_file\n";
    foreach my $gc1 (@gc) {
        my $tot = 0;
        foreach my $gc2 (@gc) {
            $tot += $gct{$gc1}{$gc2} if (defined $gct{$gc1}{$gc2});
        }
        if ($tot > 0) {
            foreach my $gc2 (@gc) {
                my $cnt = 0;
                $cnt = $gct{$gc1}{$gc2} if (defined $gct{$gc1}{$gc2});
                my $frq = sprintf ("%.8f", $cnt / $tot);
	            print G "$gc1\t$gc2\t$frq\t$cnt\n";
	        }
	    }
	    else {
	        foreach my $gc2 (@gc) {
                my $cnt = 0;
                my $frq = sprintf ("%.8f", 1 / @gc);
	            print G "$gc1\t$gc2\t$frq\t$cnt\n";
	        }
	    }
    }
    close G;
    %kmer = ();
    %gct  = ();
}

sub createKmer {
    # fill a vector with all k-mer combinations
    my $k = shift @_; $k--;
	my @old = @_;
	my @new = ();
	if ($k < 1) {
		return @old;
	}
	else {
		foreach my $e (@old) {
			foreach my $n (@dna) {
				push @new, "$e$n"; # add new element
			}
		}
		createKmer($k, @new); # recursive!
	}
}

sub profileRepeats {
    # call the repeat parsers and write the final table
    warn "profiling repeats\n" if (defined $verbose);
    profileTRF() if (defined $trf);
    profileRM()  if (defined $repeat);
    
    my $file = "$model.repeats.W$win.data";
    warn "writing repeats info in \"$file\"\n";
    open R, ">$file" or die "cannot open $file\n";
    foreach my $gc (@gc) {
        print R "#GC=$gc\n";
        foreach my $rep (@{ $repeat{$gc} }) {
            print R "$rep\n";
        }
    }
    close R;
}

sub profileTRF {
    # parse TRF output, remove overlapping repeats 
    warn "  parsing TRF files\n" if (defined $verbose);
    foreach my $file (@trf) {
        if (defined $exclude) {
            next if ($file =~ m/$exclude/);
        }
        my $last_ini = -1;
        my $last_end = -1;
        my $last_div = -1;
        my $last_gc  = -1;
        open T, "$file" or die "cannot open $file\n";
        while (<T>) {
            chomp;
            my @line      = split (/\t/, $_);
            my $seq_id    = $line[0];
            my $ini       = $line[1];
            my $end       = $line[2];
            my $period    = $line[5];
            my $div       = 100 - $line[7];
            my $indel     = $line[8];
            my $consensus = $line[-1];
            my $label     = "SIMPLE:$consensus:$period:$div:$indel";
            next unless (defined $seq{$seq_id});
            next if (checkGene($seq_id, $ini, $end));
            
            # Check for overlaping repeats
            if ($ini >= $last_ini and $ini <= $last_end) {
                if ($div >= $last_div) {
                    next; # because last repeat has less variation
                }
                else {
                    pop @{ $repeat{$last_gc} }; # last repeat removal
                }
            }
             
            my $gc    = getBinGC($seq_id, $ini);
            push @{ $repeat{$gc} }, $label;
            $last_ini = $ini;
            $last_end = $end;
            $last_div = $div;
            $last_gc  = $gc;
        }
        close T;
    }
}

sub profileRM {
    # parse TRF output, mix spliced repeats
    warn "  parsing RepeatMasker files\n" if (defined $verbose);
    my %repdata = ();
    foreach my $file (@repeat) {
        if (defined $exclude) {
            next if ($file =~ m/$exclude/);
        }
        open T, "$file" or die "cannot open $file\n";
        while (<T>) {
            chomp;
            s/^\s+//;
            next unless (m/^\d+/);
            next if (m/Simple_repeat|Low_complexity|Unknown|Satellite/);
            my @line      = split (/\s+/, $_);
            my $seq_id    = $line[4];
            my $ini       = $line[5];
            my $end       = $line[6];
            my $div       = $line[1];
            my $ins       = $line[2];
            my $del       = $line[3];
            my $dir       = $line[8];
            my $type      = $line[9];
            my $fam       = $line[10];
            my $rini      = $line[11];
            my $rend      = $line[12];
            my $rid       = "$seq_id:$fam:$line[-1]";
            next if ($type =~ m/\(/);
            if ($dir eq 'C') {
                $rini = $line[13];
                $rend = $line[12];
                $dir  = '-';
            }
            my $label     = "$type:$fam:$dir:$div:$ins:$del:$rini:$rend";
            next unless (defined $seq{$seq_id});
            next if (checkGene($seq_id, $ini, $end));
                
            if (defined $repdata{$rid}) {
                $repdata{$rid}{'label'} .= ";$div:$ins:$del:$rini:$rend";
                $repdata{$rid}{'end'}    = $end;

            }
            else {
                $repdata{$rid}{'label'}  = $label;
                $repdata{$rid}{'ini'}    = $ini;
                $repdata{$rid}{'end'}    = $end;
                $repdata{$rid}{'seq_id'} = $seq_id;
            }
        }
        close T;
    }
    
    foreach my $rid (keys %repdata) {
        my $gc = getBinGC($repdata{$rid}{'seq_id'}, int(($repdata{$rid}{'end'} - $repdata{$rid}{'ini'})/2));
        push @{ $repeat{$gc} }, $repdata{$rid}{'label'};
    }
    %repdata = ();
}

sub checkGene {
    my ($c, $i, $e) = @_;
    my $res = undef;
    foreach my $coord (@{ $genes{$c} }) {
        my ($gi, $ge) = split (/-/, $coord);
        if ($i >= $gi and $e <= $ge) {
            $res = 1;
            last;
        }
        last if ($gi > $e);
    }
    return $res;
}

sub writeMaskSeq {
    # write the masked sequence: [ACGT]=effective bases, N=ambiguous bases,
    # S=simple repeat bases, R=interspersed repeats, X=gene bases
    my $file = shift @_;
    warn "writing sequence in $file\n" if (defined $verbose);
    my $good_bases = 0;
    my $null_bases = 0;
    my $rep_bases  = 0;
    my $low_bases  = 0;
    my $func_bases = 0;
    my $tot_bases  = 0;
    open F, ">$file" or die "cannot write $file\n";
    foreach my $id (keys %seq) {
        my $seq = $seq{$id};
        print F ">$id\n";
        while ($seq) {
            my $s = substr ($seq, 0, 70);
            print F "$s\n";
            substr ($seq, 0, 70) = '';
            $good_bases += $s =~ tr/ACGT//;
            $null_bases += $s =~ tr/N//;
            $rep_bases  += $s =~ tr/R//;
            $low_bases  += $s =~ tr/S//;
            $func_bases += $s =~ tr/X//;
            $tot_bases  += length $s;
        }
    }
    close F;
	warn "compressing $file\n" if (defined $verbose);
	system ("gzip -f $file");
    warn <<__RES__
    Effective bases            = $good_bases
    Null bases (N)             = $null_bases
    Interspersed repeats bases = $rep_bases
    Low complexity bases       = $low_bases
    Functional bases (genes)   = $func_bases
    Total bases                = $tot_bases
__RES__
if (defined $verbose);
}

sub writeModelInfo {
    # write basic model description and related files created
    my $tot_bases  = 0;
    my $tot_null   = 0;
    my $tot_repeat = 0;
    my $tot_simple = 0;
    my $tot_good   = 0;
    # counting bases
    foreach my $id (keys %seq) {
        $tot_bases += length $seq{$id};
        $tot_null  += $seq{$id} =~ tr/N//;
        $tot_good  += $seq{$id} =~ tr/ACGT//;
    }
    foreach my $gc (@gc) {
        foreach my $rep (@{ $repeat{$gc} }) {
            if ($rep =~ m/^SIMPLE/) { $tot_simple++; }
            else                    { $tot_repeat++; }
        }
    }
    
    open  M, ">$model.model" or die "cannot open $model.model\n";
    print M "model=$model\n";
    print M "bases=$tot_bases\n";
    print M "intergenic=$tot_good\n";
    print M "undefined=$tot_null\n";
    print M "num_repeat=$tot_repeat\n"                  unless (defined $no_repeat_table);
    print M "num_simple=$tot_simple\n"                  unless (defined $no_repeat_table);
    close M;   
}

sub removeTmp {
    warn "removing temporary files\n" if (defined $verbose);
    my %keep = ();
    if (defined $keep_dw_files) {
        $keep{ $files{$model}{'FAS'} } = 1;
        $keep{ $files{$model}{'RMO'} } = 1;
        $keep{ $files{$model}{'TRF'} } = 1;
        $keep{ $files{$model}{'GEN'} } = 1;
    }
    
    my @dirs = qw/fasta RM TRF/;
    
    unlink ($files{$model}{'GEN'}) unless (defined $keep{ $files{$model}{'GEN'} });
    foreach my $dir (@dirs) {
        chdir $dir;
        opendir D, "." or die "cannot open directory $dir\n";
        while (my $file = readdir D) {
            next if (defined $keep{$file});
            next if ($file eq '.' or $file eq '..'); 
            system ("rm -rf $file");
        }
        closedir D;
        chdir "..";
    }
}

sub calcGC {
    my $seq = shift @_;
    $seq =~ s/[^ACGTacgt]//;
    my $len = length $seq;
    return 'NA' if ($len < 1);
    my $ngc = $seq =~ tr/CGcg/CGcg/;
    my $gc  = $ngc / $len;
    return $gc;
}

sub classGC {
    # GC ranges must be equal to @gc values or weird stuff will happen
    my $gc = shift @_;
    return $gc if ($gc eq 'NA');
    
    if    ($gc <= 0.05) { $gc =  '0-5'  ; }
    elsif ($gc <= 0.10) { $gc =  '5-10' ; }   
    elsif ($gc <= 0.15) { $gc = '10-15' ; }
    elsif ($gc <= 0.20) { $gc = '15-20' ; }
    elsif ($gc <= 0.25) { $gc = '20-25' ; }
    elsif ($gc <= 0.30) { $gc = '25-30' ; }
    elsif ($gc <= 0.35) { $gc = '30-35' ; }
    elsif ($gc <= 0.40) { $gc = '35-40' ; }
    elsif ($gc <= 0.45) { $gc = '40-45' ; }
    elsif ($gc <= 0.50) { $gc = '45-50' ; }
    elsif ($gc <= 0.55) { $gc = '50-55' ; }
    elsif ($gc <= 0.60) { $gc = '55-60' ; }
    elsif ($gc <= 0.65) { $gc = '60-65' ; }
    elsif ($gc <= 0.70) { $gc = '65-70' ; }
    elsif ($gc <= 0.75) { $gc = '70-75' ; }    
    elsif ($gc <= 0.80) { $gc = '75-80' ; }
    elsif ($gc <= 0.85) { $gc = '80-85' ; }
    elsif ($gc <= 0.90) { $gc = '85-90' ; }
    elsif ($gc <= 0.95) { $gc = '90-95' ; }
    else                { $gc = '95-100'; }
    
    return $gc;
}

sub calcBinGC {
    warn "computing GC bins\n" if (defined $verbose);
	while ( ($seq_id, $seq) = each %seq) {
		my $len  = length $seq;
		my $gc   = undef;
		my $last = $len - $win;
		my @gc   = ();
		for (my $i = 0; $i <= $last; $i += $win) {
			my $s  = substr ($seq, $i, $win);
			$gc = calcGC($s);
			push @gc, $gc;
		}
		push @gc, $gc; # last fragment with length < win
		
		my $blk = int(($binsize / $win) / 2);
		for (my $i = 0; $i <= $#gc; $i++) {
		    my $sum = 0;
		    my $num = 0;
		       $gc  = 'NA';
		    my $ini = $i - $blk; 
		       $ini = 0    if ($ini <    0);
		    my $end = $i + $blk; 
		       $end = $#gc if ($end > $#gc);
		    for (my $j = $ini; $j <= $end; $j++) {
		        next if ($gc[$j] eq 'NA');
		        $sum += $gc[$j];
		        $num++;
		    }
		    if ($num > 0) {
		        $gc = classGC($sum / $num); # average GC
		    }
		    push @{ $bingc{$seq_id} }, $gc;
		}
	}			
}

sub getBinGC {
	my $id  = shift @_;
	my $pos = shift @_;
	my $gc  = $bingc{$id}[int($pos / $win)];
	return $gc;
}

sub checkRevComp {
    my $w =  shift @_;
    my $r =  reverse $w;
       $r =~ tr/ACGTacgt/TGCAtgca/;
    my @w = ($w, $r);
       @w = sort @w;
       $w = shift @w;
    return $w;
}

sub loadFiles {
    # Human
    $files{'hg19'   }{'FAS'} = 'chromFa.tar.gz';
    $files{'hg19'   }{'RMO'} = 'chromOut.tar.gz';
    $files{'hg19'   }{'TRF'} = 'chromTrf.tar.gz';
    $files{'hg19'   }{'GEN'} = 'ensGene.txt.gz';
    
    $files{'hg18'   }{'FAS'} = 'chromFa.tar.gz';
    $files{'hg18'   }{'RMO'} = 'chromOut.tar.gz';
    $files{'hg18'   }{'TRF'} = 'chromTrf.tar.gz';
    $files{'hg18'   }{'GEN'} = 'ensGene.txt.gz';
    
    # Cat
    $files{'felCat4'}{'FAS'} = 'felCat4.fa.gz';
    $files{'felCat4'}{'RMO'} = 'felCat4.fa.out.gz';
    $files{'felCat4'}{'TRF'} = 'felCat4.trf.bed.gz';
    $files{'felCat4'}{'GEN'} = 'refGene.txt.gz';
    
    # Chicken
    $files{'galGal3'}{'FAS'} = 'chromFa.tar.gz';
    $files{'galGal3'}{'RMO'} = 'chromOut.tar.gz';
    $files{'galGal3'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'galGal3'}{'GEN'} = 'ensGene.txt.gz';
    
    # Chimpanzee
    $files{'panTro3'}{'FAS'} = 'panTro3.fa.gz';
    $files{'panTro3'}{'RMO'} = 'panTro3.fa.out.gz';
    $files{'panTro3'}{'TRF'} = 'panTro3.trf.bed.gz';
    $files{'panTro3'}{'GEN'} = 'refGene.txt.gz';
    
    # Cow
    $files{'bosTau4'}{'FAS'} = 'bosTau4.fa.gz';
    $files{'bosTau4'}{'RMO'} = 'bosTau4.fa.out.gz';
    $files{'bosTau4'}{'TRF'} = 'bosTau4.trf.bed.gz';
    $files{'bosTau4'}{'GEN'} = 'ensGene.txt.gz';
    
    # Dog
    $files{'canFam2'}{'FAS'} = 'chromFa.tar.gz';
    $files{'canFam2'}{'RMO'} = 'chromOut.tar.gz';
    $files{'canFam2'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'canFam2'}{'GEN'} = 'ensGene.txt.gz';
    
    # Elephant
    $files{'loxAfr3'}{'FAS'} = 'loxAfr3.fa.gz';
    $files{'loxAfr3'}{'RMO'} = 'loxAfr3.fa.out.gz';
    $files{'loxAfr3'}{'TRF'} = 'loxAfr3.trf.bed.gz';
    $files{'loxAfr3'}{'GEN'} = 'ensGene.txt.gz';
    
    # Fugu
    $files{'fr2'    }{'FAS'} = 'chromFa.tar.gz';
    $files{'fr2'    }{'RMO'} = 'chromOut.tar.gz';
    $files{'fr2'    }{'TRF'} = 'chromTrf.tar.gz';
    $files{'fr2'    }{'GEN'} = 'ensGene.txt.gz'; 
    
    # Guinea Pig
    $files{'cavPor3'}{'FAS'} = 'cavPor3.fa.gz';
    $files{'cavPor3'}{'RMO'} = 'cavPor3.fa.out.gz';
    $files{'cavPor3'}{'TRF'} = 'cavPor3.trf.bed.gz';
    $files{'cavPor3'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Horse
    $files{'equCab2'}{'FAS'} = 'chromFa.tar.gz';
    $files{'equCab2'}{'RMO'} = 'chromOut.tar.gz';
    $files{'equCab2'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'equCab2'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Lizard
    $files{'anoCar2'}{'FAS'} = 'anoCar2.fa.gz';
    $files{'anoCar2'}{'RMO'} = 'anoCar2.fa.out.gz';
    $files{'anoCar2'}{'TRF'} = 'anoCar2.trf.bed.gz';
    $files{'anoCar2'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Marmoset
    $files{'calJac3'}{'FAS'} = 'calJac3.fa.gz';
    $files{'calJac3'}{'RMO'} = 'calJac3.fa.out.gz';
    $files{'calJac3'}{'TRF'} = 'calJac3.trf.bed.gz';
    $files{'calJac3'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Medaka
    $files{'oryLat2'}{'FAS'} = 'oryLat2.fa.gz';
    $files{'oryLat2'}{'RMO'} = 'oryLat2.fa.out.gz';
    $files{'oryLat2'}{'TRF'} = 'oryLat2.trf.bed.gz';
    $files{'oryLat2'}{'GEN'} = 'ensGene.txt.gz'; 
        
    # Mouse
    $files{'mm9'    }{'FAS'} = 'chromFa.tar.gz';
    $files{'mm9'    }{'RMO'} = 'chromOut.tar.gz';
    $files{'mm9'    }{'TRF'} = 'chromTrf.tar.gz';
    $files{'mm9'    }{'GEN'} = 'ensGene.txt.gz'; 

    # Oposum
    $files{'monDom5'}{'FAS'} = 'chromFa.tar.gz';
    $files{'monDom5'}{'RMO'} = 'chromOut.tar.gz';
    $files{'monDom5'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'monDom5'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Orangutan
    $files{'ponAbe2'}{'FAS'} = 'chromFa.tar.gz';
    $files{'ponAbe2'}{'RMO'} = 'chromOut.tar.gz';
    $files{'ponAbe2'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'ponAbe2'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Panda
    $files{'ailMel1'}{'FAS'} = 'ailMel1.fa.gz';
    $files{'ailMel1'}{'RMO'} = 'ailMel1.fa.out.gz';
    $files{'ailMel1'}{'TRF'} = 'ailMel1.trf.bed.gz';
    $files{'ailMel1'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Pig
    $files{'susScr2'}{'FAS'} = 'chromFa.tar.gz';
    $files{'susScr2'}{'RMO'} = 'chromOut.tar.gz';
    $files{'susScr2'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'susScr2'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Platypus
    $files{'ornAna1'}{'FAS'} = 'ornAna1.fa.gz';
    $files{'ornAna1'}{'RMO'} = 'ornAna1.fa.out.gz';
    $files{'ornAna1'}{'TRF'} = 'ornAna1.trf.bed.gz';
    $files{'ornAna1'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Rabbit
    $files{'oryCun2'}{'FAS'} = 'oryCun2.fa.gz';
    $files{'oryCun2'}{'RMO'} = 'oryCun2.fa.out.gz';
    $files{'oryCun2'}{'TRF'} = 'oryCun2.trf.bed.gz';
    $files{'oryCun2'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Rat
    $files{'rn4'    }{'FAS'} = 'chromFa.tar.gz';
    $files{'rn4'    }{'RMO'} = 'chromOut.tar.gz';
    $files{'rn4'    }{'TRF'} = 'chromTrf.tar.gz';
    $files{'rn4'    }{'GEN'} = 'ensGene.txt.gz'; 
    
    # Rhesus
    $files{'rheMac2'}{'FAS'} = 'chromFa.tar.gz';
    $files{'rheMac2'}{'RMO'} = 'chromOut.tar.gz';
    $files{'rheMac2'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'rheMac2'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Sheep
    $files{'oviAri1'}{'FAS'} = 'oviAri1.fa.gz';
    $files{'oviAri1'}{'RMO'} = 'oviAri1.fa.out.gz';
    $files{'oviAri1'}{'TRF'} = 'oviAri1.trf.bed.gz';
    $files{'oviAri1'}{'GEN'} = 'refGene.txt.gz'; 
    
    # Stickleback
    $files{'gasAcu1'}{'FAS'} = 'chromFa.tar.gz';
    $files{'gasAcu1'}{'RMO'} = 'chromOut.tar.gz';
    $files{'gasAcu1'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'gasAcu1'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Tetraodon
    $files{'tetNig2'}{'FAS'} = 'chromFa.tar.gz';
    $files{'tetNig2'}{'RMO'} = 'chromOut.tar.gz';
    $files{'tetNig2'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'tetNig2'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # X. tropicalis
    $files{'xenTro2'}{'FAS'} = 'xenTro2.fa.gz';
    $files{'xenTro2'}{'RMO'} = 'xenTro2.rmsk.out.gz';
    $files{'xenTro2'}{'TRF'} = 'xenTro2.trf.bed.gz';
    $files{'xenTro2'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Zebra finch
    $files{'taeGut1'}{'FAS'} = 'chromFa.tar.gz';
    $files{'taeGut1'}{'RMO'} = 'chromOut.tar.gz';
    $files{'taeGut1'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'taeGut1'}{'GEN'} = 'ensGene.txt.gz'; 
    
    # Zefrafish
    $files{'danRer7'}{'FAS'} = 'danRer7.fa.gz';
    $files{'danRer7'}{'RMO'} = 'danRer7.fa.out.gz';
    $files{'danRer7'}{'TRF'} = 'danRer7.trf.bed.gz';
    $files{'danRer7'}{'GEN'} = 'ensGene.txt.gz';
    
    # C. intestinalis
    $files{'ci2'    }{'FAS'} = 'ScaffoldFa.zip';
    $files{'ci2'    }{'RMO'} = 'Scaffold.out.zip';
    $files{'ci2'    }{'TRF'} = 'ScaffoldTrf.zip';
    $files{'ci2'    }{'GEN'} = 'ensGene.txt.gz';
    
    # S. purpuratus
    $files{'strPur2'}{'FAS'} = 'strPur2.fa.gz';
    $files{'strPur2'}{'RMO'} = 'strPur2.fa.out.gz';
    $files{'strPur2'}{'TRF'} = 'strPur2.trf.bed.gz';
    $files{'strPur2'}{'GEN'} = 'refGene.txt.gz';
    
    # A. gambiae
    $files{'anoGam1'}{'FAS'} = 'chromFa.tar.gz';
    $files{'anoGam1'}{'RMO'} = 'chromOut.tar.gz';
    $files{'anoGam1'}{'TRF'} = 'chromTrf.tar.gz';
    $files{'anoGam1'}{'GEN'} = 'ensGene.txt.gz';
    
    # A. mellifera
    $files{'apiMel2'}{'FAS'} = 'GroupFa.zip';
    $files{'apiMel2'}{'RMO'} = 'GroupOut.zip';
    $files{'apiMel2'}{'TRF'} = 'GroupTrf.zip';
    $files{'apiMel2'}{'GEN'} = 'ensGene.txt.gz';
    
    # D. melanogaster
    $files{'dm3'    }{'FAS'} = 'chromFa.tar.gz';
    $files{'dm3'    }{'RMO'} = 'chromOut.tar.gz';
    $files{'dm3'    }{'TRF'} = 'chromTrf.tar.gz';
    $files{'dm3'    }{'GEN'} = 'ensGene.txt.gz';
    
    # C. elegans
    $files{'ce3'    }{'FAS'} = 'chromFa.tar.gz';
    $files{'ce3'    }{'RMO'} = 'chromOut.tar.gz';
    $files{'ce3'    }{'TRF'} = 'chromTrf.tar.gz';
    $files{'ce3'    }{'GEN'} = 'ensGene.txt.gz';
}
