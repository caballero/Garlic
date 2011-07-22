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

The model names are according to UCSC:
  - hg19    -> human genome release 19 (GRCh37)
  - hg18    -> human genome release 18 (GRCh36)
  - mm9     -> mouse genome release 9
  - equCab2 -> horse genome release 2
  - ...

This script will download the required files and process it. Caution, the files
are really big and the download time could be large.

Also you can use your own sequences and annotations to create a model:

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
my $write_mask_seq  = undef;      # Write fasta flag

# Fetch options
GetOptions(
    'h|help'             => \$help,
    'v|verbose'          => \$verbose,
    'm|model=s'          => \$model,
    'd|dir:s'            => \$dir,
    'k|kmer:i'           => \$kmer,
    'w|win:i'            => \$win,
    'f|fasta:s'          => \$fasta,
    'r|repeat:s'         => \$repeat,
    't|trf:s'            => \$trf,
    'g|gene:s'           => \$gene,
    'e|exclude:s'        => \$exclude,
    'rm_tmp'             => \$rm_tmp,
    'mask_repeat|mR'     => \$mask_repeat,
    'mask_trf|mT'        => \$mask_trf,
    'no_mask_gene|nG'    => \$no_mask_gene,
    'no_repeat_table|nR' => \$no_repeat_table,
    'no_kmer_table|nK'   => \$no_kmer_table,
    'write_mask_seq'     => \$write_mask_seq
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
my %simple   = ();
my %repeat   = ();
my @dna      = qw/A C G T/;
my ($seq, $ss, $seq_id, $ini, $end, $len);
my @gc = qw/0-10 10-20 20-30 30-40 40-50 50-60 60-70 70-80 80-90 90-100/;

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
getUCSC_repeat() unless (defined $repeat or defined $no_repeat_table);
getUCSC_trf()    unless (defined $trf    or defined $no_repeat_table);
getUCSC_gene()   unless (defined $gene   or defined $no_mask_gene);

# Creating list of files to process
@fasta   = split (/,/, $fasta);
@repeat  = split (/,/, $repeat);
@trf     = split (/,/, $trf);
@gene    = split (/,/, $gene);
$exclude =~ s/,/|/g if (defined $exclude);

# Load sequences and mask them
readFasta();
maskGene() unless (defined $no_mask_gene);
maskRepeat()   if (defined $mask_repeat);
maskTRF()      if (defined $mask_trf);
writeMaskSeq("$model.masked.fa") if (defined $write_mask_seq);

# Create the K-mer/Window and GCt tables
profileSeqs() unless (defined $no_kmer_table);

# Create the Repeat Table
profileRepeats() unless (defined $no_repeat_table);

# Clean up
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
    $gene = 'ensGene.txt.gz';
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
    my $target_file = 'chromTrf.tar.gz';
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
    
    warn "unpacking TAR\n" if (defined $verbose);
    system ("$unpack $target_file");
    chdir '..';
    warn "searching TRF files\n" if (defined $verbose);
    $trf = searchFiles('.bed$', 'TRF');
    die "cannot find TRF files!" unless (defined $trf);
}

sub getUCSC_repeat {
    my $target_file = 'chromOut.tar.gz';
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
    
    warn "unpacking TAR\n" if (defined $verbose);
    system ("$unpack $target_file");
    chdir '..';
    warn "searching RM files\n" if (defined $verbose);
    $repeat = searchFiles('.out$', 'RM');
    die "cannot find RM files!" unless (defined $repeat);
}

sub getUCSC_fasta {
    my $target_file = 'chromFa.tar.gz';
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
    
    warn "unpacking TAR\n" if (defined $verbose);
    system ("$unpack $target_file");
    chdir '..';
    warn "searching fasta files\n" if (defined $verbose);
    $fasta = searchFiles('.fa$', 'fasta');
    die "cannot find fasta files!" unless (defined $fasta);
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
    warn "masking repeats\n" if (defined $verbose);
    foreach my $file (@repeat) {
        if (defined $exclude) {
            next if ($file =~ m/$exclude/);
        }
        my $fileh = $file;
        $fileh = "gunzip  -c $file | " if ($file =~ m/\.gz$/);
        $fileh = "bunzip2 -c $file | " if ($file =~ m/\.bz2$/);
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
    warn "masking simple repeats\n" if (defined $verbose);
    foreach my $file (@trf) {
        if (defined $exclude) {
            next if ($file =~ m/$exclude/);
        }
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
            substr ($seq{$seq_id}, $ini, $len) = 'S' x $len;
        }
        close FH;
    }
}

sub maskGene {
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
            $len = $end - $ini - 1;
            substr ($seq{$seq_id}, $ini, $len) = 'X' x $len;
        }
        close FH;
    }
}

sub profileSeqs {
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
	    for (my $i = 0; $i <= (length $seq) - $win; $i += $win) {
	        $ss = substr ($seq, $i, $win);
	        my $num_N = $ss =~ tr/N/N/;
	        my $frq_N = $num_N / (length $ss);
	        next if ($frq_N > 0.3); # No more than 30% bad bases 
	        $bp_slices += (length $ss) - $num_N; # Effective bases
	        
	        # GC in the window and transition
	        my $gc = calcGC($ss);
	        $gct{$last_gc}{$gc}++ if (defined $last_gc);
	        $last_gc = $gc;
	       
	        # Kmer count
	        for (my $j = 0; $j <= $win - $kmer; $j++) {
	            my $kmer = substr ($ss, $j, $kmer);
	            $kmer{$gc}{$kmer}++;
	        }
        }
    }
    
    warn "  $bp_slices bases analyzed\n" if (defined $verbose);
    my $kmer_file = "$model.kmer.K$kmer.W$win.data";
    warn "writing k-mer profile in \"$kmer_file\"\n" if (defined $verbose);
    open K, ">$kmer_file" or die "cannot write file $kmer_file\n";
    foreach my $gc (@gc) {
	    print K "#GC=$gc\n";
	    foreach my $word (@kmer) {
	        my $tot = 0;
	        foreach my $b (@dna) {
	            $tot += $kmer{$gc}{"$word$b"} if (defined $kmer{$gc}{"$word$b"});
	        }
	        if ($tot > 0) {
	            foreach my $b (@dna) {
	                my $cnt = 0;
	                $cnt = $kmer{$gc}{"$word$b"} if (defined $kmer{$gc}{"$word$b"});
	                my $frq = sprintf ("%.8f", $cnt / $tot);
	                print K "$word$b\t$frq\t$cnt\n";
	            }
	        }
	        else {
	            foreach my $b (@dna) {
	                my $frq = sprintf ("%.8f", 1 / @dna);
	                print K "$word$b\t$frq\t0\n";
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
		createKmer($k, @new);
	}
}

sub profileRepeats {
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
    %repeat = ();
}

sub profileTRF {
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
            
            # Check for overlaping repeats
            if ($ini >= $last_ini and $ini <= $last_end) {
                if ($div >= $last_div) {
                    next; # because last repeat has less variation
                }
                else {
                    pop @{ $repeat{$last_gc} }; # last repeat removal
                }
            }
             
            my $left      = substr ($seq{$seq_id}, $ini - $win, $win);
            my $right     = substr ($seq{$seq_id}, $end, $win);
            my $gc        = calcGC("$left$right");
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
            next if (m/Simple_repeat|Low_complexity|Unknown/);
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
            if ($dir eq 'C') {
                $rini = $line[13];
                $rend = $line[12];
                $dir  = '-';
            }
            my $label     = "$type:$fam:$dir:$div:$ins:$del:$rini:$rend";
            next unless (defined $seq{$seq_id});
            my $left      = substr ($seq{$seq_id}, $ini - $win, $win);
            my $right     = substr ($seq{$seq_id}, $end, $win);
                
            if (defined $repdata{$rid}) {
                $repdata{$rid}{'rseq'}   = $right;
                $repdata{$rid}{'label'} .= ";$div:$ins:$del:$rini:$rend";
            }
            else {
                $repdata{$rid}{'label'} = $label;
                $repdata{$rid}{'lseq'}  = $left;
                $repdata{$rid}{'rseq'}  = $right;
            }
        }
        close T;
    }
    
    foreach my $rid (keys %repdata) {
        my $gc = calcGC($repdata{$rid}{'lseq'} . $repdata{$rid}{'rseq'});
        push @{ $repeat{$gc} }, $repdata{$rid}{'label'};
    }
}

sub writeMaskSeq {
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
    # GC ranges must be equal to @gc values or weird stuff will happen
    my $seq = shift @_;
    $seq =~ s/[^ACGTacgt]//;
    my $len = length $seq;
    return 'NA' if ($len < 1);
    my $ngc = $seq =~ tr/CGcg/CGcg/;

    my $gc  = $ngc / $len;
    if    ($gc <= 0.10) { $gc =  '0-10' ; }
    elsif ($gc <= 0.20) { $gc = '10-20' ; }   
    elsif ($gc <= 0.30) { $gc = '20-30' ; }
    elsif ($gc <= 0.40) { $gc = '30-40' ; }
    elsif ($gc <= 0.50) { $gc = '40-50' ; }
    elsif ($gc <= 0.60) { $gc = '50-60' ; }   
    elsif ($gc <= 0.70) { $gc = '60-70' ; }
    elsif ($gc <= 0.80) { $gc = '70-80' ; }
    elsif ($gc <= 0.90) { $gc = '80-90' ; }
    else                { $gc = '90-100'; }
    
    return $gc;
}
    
