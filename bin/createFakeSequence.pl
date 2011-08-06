#!/usr/bin/perl

=head1 NAME 

createFakeSequence.pl

=head1 DESCRIPTION

Perl script to create a sequence similar to the intergenic regions from
the model used. The models are described in parts:

a) The composition background with Markov models from 0 to 5 in fixed windows.
b) The repeats parsed from the alignements described by RepeatMasker and 
consensus bases acording to RepBase.
c) Transitional frecuencies for fixed windows in GC% classes (0-90).

To create the new sequence first we need to retrieve the elements to insert 
and calculate the length of the base DNA. Then, the elements are "bombarded"
in random positions.

The main idea is to recreate a sequence generator without limits of size and
sequences as real as the intergenic regions of a genome.

=head1 USAGE

perl createFakeSequence.pl [OPTIONS]

=head1 EXAMPLES

** Examples

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
use List::Util qw/shuffle/;

## Global variables
my $seq     =       ''; # The sequence itself
my $model   =    undef; # Model to use
my $kmer    =        4; # HMM size to use
my $nrep    =    undef; # Number of repeats to insert
my $nsim    =    undef; # Number of simple repeats to insert
my $out     =   'fake'; # Filename to use
my $size    =    undef; # Final size of the sequence
my $mask    =    undef; # Flag to mask repeats
my $win     =     1000; # Window size for sequence GC transition
my $help    =    undef; # Variable to activate help
my $debug   =    undef; # Variable to activate verbose mode
my %model   =       (); # Hash to store model parameters
my %simple  =       (); # Hash to store repeats info
my %repeat  =       (); # Hash to store repeats info
my @inserts =       (); # Array to store the insert sequences data
my %gct     =       (); # Hash for GC transitions GC(n-1) -> GC(n)
my %elemk   =       (); # Hash for kmers probabilities
my %gc      =       (); # Hash for GC content probabilities
my @classgc =       (); # Array for class GC
my %classgc =       (); # Hash for class GC
my $mingc   =       10; # Minimal GC content to use
my $maxgc   =      100; # Maximal GC content to use
my $dir     = './data'; # Path to models/RebBase directories
my %rep_seq =       (); # Hash with consensus sequences from RepBase
my $max_cyc =     1000; # Max number of cycles in loops
my $mut_cyc =       10;
my @dna     = qw/A C G T/; # yes, the DNA alphabet

## Parameters extraction
usage() if (!GetOptions( 
				'help|h'    	=> \$help,
				'model|o=s' 	=> \$model,
				'length|l=s'  	=> \$size,
				'name|n=s'   	=> \$out,
				'kmer|k:i'  	=> \$kmer,
				'win|w:i'   	=> \$win,
				'repeats|r:i'  	=> \$nrep,
				'simple|s:i'  	=> \$nsim,
				'mingc|g:i'     => \$mingc,
				'maxgc|c:i'     => \$maxgc,
				'mask|m'     	=> \$mask,
				'verbose|v'   	=> \$debug,
				'dir:s'         => \$dir
				)
);

usage() if (defined $help);
usage() unless (defined $model and defined $size);

# Loading model parameters
readConfig("$dir/$model/$model.model");
$model{'gct_file'}    = "$model.GCt.W$win.data";
$model{'kmer_file'}   = "$model.kmer.K$kmer.W$win.data";
$model{'repeat_file'} = "$model.repeats.W$win.data";
$model{'repbase'}     = "$dir/RepBase/RepBase16.06.fa.gz";    # point to RepBase fasta file

# GC classes creation
for (my $i = $mingc; $i <= $maxgc; $i += 5) { 
	push @classgc, $i;
	$classgc{$i}++;
}

print "Generating a $size sequence with $model model, output in $out\n" if (defined $debug);

# Checking the size (conversion of symbols)
$size = checkSize($size);
errorExit("$size isn't a number") unless ($size > 1);

# Loading background models
loadGCt("$dir/$model/" . $model{'gct_file'});
print "GC transitions in ", $model{'gct_file'}, " loaded\n" if(defined $debug);

loadKmers("$dir/$model/" . $model{'kmer_file'});
print "k-mers in ", $model{'kmer_file'}, " loaded\n" if(defined $debug);

# Number of simple repeats to use
unless (defined $nsim) {
	$nsim = calcInsertNum($size, $model{'num_simple'} / $model{'intergenic'});
}
print "$nsim simple repeats to select\n" if(defined $debug);

# Number of interspearsed repeats to use
unless (defined $nrep) {
	$nrep = calcInsertNum($size, $model{'num_repeat'} / $model{'intergenic'});
}
print "$nrep interspersed repeats to select\n" if(defined $debug);

if ($nrep > 0) {
	loadRepeatConsensus($model{'repbase'});
}

if (($nrep + $nsim) > 0) {
    loadRepeats("$dir/$model/" . $model{'repeat_file'});
}

# Generation of base sequence
my $fgc    = newGC();
my @fseeds = keys %{ $elemk{$fgc} };
my $fseed  = $fseeds[int(rand @fseeds)];

$seq = createSeq($kmer, $fgc, $size, $win, $fseed);
print "Base sequence generated ($size bases)\n" if(defined $debug);

# Adding new elements
if ($nsim + $nrep > 0) {
	$seq = insertElements($seq);
	print "Inserted [$nsim + $nrep] elements in base sequence\n" if(defined $debug);
}

print "Verifying sequence final size\n" if(defined $debug);
$seq = checkSeqSize($size, $seq);

print "Generated a sequence with ", length $seq, " bases\n" if(defined $debug);

# Printing output
print "Printing outputs\n" if(defined $debug);
my $fseq = formatFasta($seq);
open  FAS, ">$out.fasta" or errorExit("cannot open $out.fasta");
print FAS  ">artificial_sequence MODEL=$model KMER=$kmer WIN=$win LENGTH=$size\n$fseq";
close FAS;

if ($nsim + $nrep > 0) {
    open  INS, ">$out.inserts" or errorExit("cannot open $out.inserts");
    print INS  "POS\tREPEAT\n";
    foreach my $rep (@inserts) {
        print INS "$rep\n";	
	}
    close INS;
}

#################################################
##      S  U  B  R  O  U  T  I  N  E  S        ##
#################################################

sub usage {
print <<__HELP__
Usage: perl intergenic.pl [--help|-h] -o MODEL -l SIZE -n OUFILE [PARAMETERS] 
Parameters:
  -o --model     Model to use (like hg19, mm9, ... etc).
  -l --length    Size in bases [kb, Mb, Gb accepted].
  -n --name      Output files to create [*.fasta and *.log].
	
Optional or automatic parameters:
  -w --win       Window size for base generation profile.     Default =  1000
  -k --kmer      Seed size to use [available: 1,2,3,4,5,6].   Default =     4
  -g --mingc     Minimal GC content to use [0,10,..,90].      Default =    10
  -c --maxgc     Maximal GC content to use [0,10,..,90].      Default =   100
  -r --repeats   Number of total repeats to insert.           Default =  Auto
  -s --simple    Number of total simple repeats to insert.    Default =  Auto
  -m --mask      Mask repeats in final sequence.              Default = False
  -v --verbose   Verbose output for debug.                    Default = False
  -h --help      Print this screen.
	
__HELP__
;
exit 1;
}

# formatFasta => break a sequence in blocks (80 col/line default)
sub formatFasta {
	my $sseq  = shift @_;
	my $col   = shift @_;
	$col ||= 70;
	my $fseq = '';
	if (length $sseq <= $col) { 
		$fseq = "$sseq\n";
	}
	else {
		while ($sseq) { 
			$fseq .= substr ($sseq, 0, $col);
			$fseq .= "\n";
			substr ($sseq, 0, $col) = '';
		}
	}
	return $fseq;
}

# defineFH => check if the file is compressed (gzip/bzip2), return the handler
sub defineFH {
    my ($fo) = @_;
    my $fh = $fo;
    $fh = "gunzip  -c $fo | " if ($fo =~ m/gz$/);
    $fh = "bunzip2 -c $fo | " if ($fo =~ m/bz2$/);
    return $fh;
}

# checkBases => verify a sequence, change unclassified bases
sub checkBases {
	my $cseq = shift @_;
	$cseq =  uc $cseq;
	$cseq =~ tr/U/T/;
	if($cseq =~ /[^ACGT]/) {
		my @cseq = split (//, $cseq);
		my @bases = ();
		for (my $i = 0; $i <= $#cseq; $i++) {
			if    ($cseq[$i] =~ /[ACGT]/) { next;                 }
			elsif ($cseq[$i] eq 'R')      { @bases = qw/A G/;     }
			elsif ($cseq[$i] eq 'Y')      { @bases = qw/C T/;     }
			elsif ($cseq[$i] eq 'M')      { @bases = qw/C A/;     }
			elsif ($cseq[$i] eq 'K')      { @bases = qw/T G/;     }
			elsif ($cseq[$i] eq 'W')      { @bases = qw/T A/;     }
			elsif ($cseq[$i] eq 'S')      { @bases = qw/C G/;     }
			elsif ($cseq[$i] eq 'B')      { @bases = qw/C T G/;   }
			elsif ($cseq[$i] eq 'D')      { @bases = qw/A T G/;   }
			elsif ($cseq[$i] eq 'H')      { @bases = qw/A T C/;   }
			elsif ($cseq[$i] eq 'V')      { @bases = qw/A C G/;   }
			else  { @bases = qw/A C G T/; } # Unknown base?
			$cseq[$i] = $bases[ int(rand(@bases)) ];
		}
		$cseq = join('', @cseq);
	}
	return $cseq;
}

# revcomp => return the reverse complementary chain
sub revcomp {
	my $rseq =  shift @_;
	$rseq    =~ tr/ACGTacgt/TGCAtgca/;
	return reverse $rseq;
}

# checkSeqSize => verify the sequence length
sub checkSeqSize {
	my $size     = shift @_;
	my $seq      = shift @_;
	my $old_size = length $seq;
	if ($old_size > $size) {
		print "Sequence too long, removing ", $old_size - $size, "bases\n" if (defined $debug);
		$seq = substr ($seq, 0, $size);
	}
	elsif ($old_size < $size) {
		my $add   = $size - $old_size;
		print "Sequence too short, adding $add bases\n" if (defined $debug);
		my $pos   = int(rand($old_size - $add));
		my $patch = substr($seq, $pos, $add);
		$seq     .= $patch;
	}
	else {
		return $seq;
	}
	checkSeqSize($size, $seq); # Recursion!
}

# calcInsertNum => return a number of inserts with some variation
sub calcInsertNum {
	my $size = shift @_;
	my $freq = shift @_;
	my $avg  = int($size * $freq);
	$avg++;
	my $num  = 0;
	my @op   = qw/plus minus none/;
	my $op   = $op[int(rand(@op))];
	my $dif  = int(rand($avg) * int(rand(2)));
	if ($op eq 'plus') { 
	    $num = $avg + $dif; 
	}
	elsif ($op eq 'minus') { 
	    $num = $avg - $dif; 
	}
	else { 
	    $num = $avg; 
	}
	$num = 0 if ($num < 0);
	return $num;
}

# readConfig => well, read the configuration file
sub readConfig {
	my $file  = shift @_;
	my $fileh = defineFH($file);
	open C, "$fileh" or errorExit("cannot open $file");
	while (<C>) {
		chomp;
		next if(m/^#/);
		next unless (m/=/);
		my ($param, $value) = split (/=/, $_);
		$model{$param} = $value;
	}
	close G;
}

# checkSize => decode kb, Mb and Gb symbols
sub checkSize{
	my $sz = shift @_;
	my $fc = 1;
	if ($sz =~ m/k/i) { 
	    $fc =  1e3;       
	}
	elsif ($sz =~ m/m/i) {
	    $fc = 1e6;
	}
	elsif ($sz =~ m/g/i) { 
	    $fc =  1e9; 
	}
    $sz  =~ s/\D//g;
    $sz *=  $fc;
	return $sz;
}

# loadGCt => load the GC transitions values
sub loadGCt {
    warn "loading GC transitions table\n" if (defined $debug);
	my $file   = shift @_;
	my $fileh  = defineFH($file);
	my %gc_sum = ();
	open G, "$fileh" or errorExit("cannot open $file");
	while (<G>) {
		chomp;
		next if (m/#/);
		my($pre, $post, $p, $num) = split (/\t/, $_);
		$num++; # give a chance to zero values
		my ($a, $b) = split (/-/, $pre);
		$pre  = $b;
		($a, $b) = split (/-/, $post);
		$post = $b;
		next unless($pre  >= $mingc and $pre  <= $maxgc);
		next unless($post >= $mingc and $post <= $maxgc);
		$gct{$pre}{$post} = $num;
		$gc_sum{$pre} += $num;
	}
	close G;
	# Adjust the probability of GC
	foreach my $pre (keys %gct) {
		foreach my $post (keys %{ $gct{$pre} }) {
			$gct{$pre}{$post} /= $gc_sum{$pre};
		}
	}
}

# loadKmers => load the kmers frequencies
sub loadKmers {
    warn "loadign kmer table\n" if (defined $debug);
	my $file  = shift @_;
	my $gc    = undef; 
	my $tot   = 0;
	my $fileh = defineFH($file);
	open K, "$fileh" or errorExit("cannot open $fileh");
	while (<K>) {
		if (m/#GC=\d+-(\d+)/) {
			$gc = $1;
		}
		else {
			chomp;
			my ($b, $f, @r) = split (/\t/, $_);
			my $v = pop @r;
			$elemk{$gc}{$b} = $f;
			$v++; # give a chance to zero values
			if (defined $classgc{$gc}) {
			    $gc{$gc} += $v;
			    $tot     += $v;
			}
		}
	}
	close K;
	# Adjust the probability of GC
	foreach $gc (keys %classgc) { 
	    $gc{$gc} /= $tot; 
	}
}

# loadRepeatConsensus => read the fasta file of RepBase consensus
sub loadRepeatConsensus {
    warn "loading interspersed repeats consensus\n" if (defined $debug);
	my $file  = shift @_;
	my $fileh = defineFH($file);
	my $rep   = '';
	open R, "$fileh" or errorExit("Cannot open $fileh");
	while (<R>) {
		chomp;
		if (m/>(\S+)/) {
			$rep = $1; 
		}
		else {
			$rep_seq{$rep} .= checkBases($_);
		}
	}
	close R;
}

# loadRepeats => read the repeats info
sub loadRepeats {
    my $file  = shift @_;
    my $fileh = defineFH($file);
    my $gc    = undef;
    my $nsim  = 0;
    my $nrep  = 0;
    open R, "$fileh" or die "cannot open $file\n";
    while (<R>) {
        chomp;
        if (m/#GC=\d+-(\d+)/) {
            $gc = $1;
        }
        else {
            next unless (defined $classgc{$gc});
            if (m/^SIMPLE/) {
                push @{ $simple{$gc} }, $_;
                $nsim++;
            }
            else {
                push @{ $repeat{$gc} }, $_;
                $nrep++;
            }
        }
    }
    close R;
    print "readed SIMPLE=$nsim REPEAT=$nrep\n"; 
}

# selPosition => find where to put an insert
sub selPosition {
	my $seq = shift @_;
	my $gc  = shift @_;
	my @pos = ();
	my %dat = ();
	for (my $i = 0; $i <= ((length $seq) - $kmer - 1); $i++) {
		my $seed = substr($seq, $i, $kmer);
		$dat{$i} = $elemk{$gc}{$seed};
	}
	foreach my $pos (sort { $dat{$a} <=> $dat{$b} } keys %dat) {
	 	push (@pos, $pos);
	}
	return @pos;
}

# addDeletions => remove bases
sub addDeletions {
	my $seq  = shift @_;
	my $ndel = shift @_;
	my $gcl  = shift @_;
	my $eval = shift @_; $eval ||= 0;
	my $tdel = 0;
	my $skip = 0;
	my @pos  = selPosition($seq, $gcl);
	
	while ($ndel > 0) {
		my $bite = 1;
		last unless(defined $pos[0]);
		my $pos  = shift @pos;
		next if ($pos < $kmer);
		next if ($pos >= length $seq);
		my $pre  = substr($seq, $pos, $kmer - 1);
		my $old  = substr($seq, $pos, $kmer);
		my $new  = substr($seq, $pos + $kmer, 1);
		if ($eval < $mut_cyc) {
			next unless(defined $elemk{$gcl}{"$pre$new"} and defined $elemk{$gcl}{"$old"});
			next if($elemk{$gcl}{"$old"} >= $elemk{$gcl}{"$pre$new"} + 1e-15);
		}
		
		substr($seq, $pos, $bite) = 'D' x $bite;
		$ndel -= $bite;
		$tdel++;
	}
	$skip = $ndel;
	print "  Added $tdel deletions ($skip skipped, GC=$gcl)\n" if(defined $debug);
	if ($skip > 0) {
		$gcl = newGC();
		$seq = addDeletions($seq, $skip, $gcl, ++$eval);
	}
	
	$seq =~ s/D//g;
	return $seq;
}

# addInsertions => add bases
sub addInsertions {
	my $seq  = shift @_;
	my $nins = shift @_;
	my $gcl  = shift @_;
	my $tins = 0;
	my @pos  = selPosition($seq, $gcl);
	#print "PreInsert: $seq\n" if(defined $debug);
	while ($nins > 0) {
		my $ins  = 1;
		last unless(defined $pos[0]);
		my $pos = shift @pos;
		next if($pos < $kmer);
		
		my $seed = substr($seq, $pos - $kmer + 1, $kmer - 1);
		my $post = substr($seq, $pos, 1);
		my $eed  = substr($seed, 1);
		my $dice = rand();
		my $n    = $dna[$#dna];
		my $p    =  0;
		unless (defined $elemk{$gcl}{"$seed$n"} and defined $elemk{$gcl}{"$eed$n$post"}) {
			print "Bad seed ($seed) in $pos ($seq)\n" if (defined $debug);
			next;
		}
		else {
			foreach my $b (@dna) {
				my $q = $p + $elemk{$gcl}{"$seed$b"};
				$n    = $b if ($dice >= $p);
				last if($dice >= $p and $dice <= ($q + 1e-15));
				$p    = $q;
			}
			next unless($dice >= $elemk{$gcl}{"$eed$n$post"} + 1e-15);
		}
		substr($seq, $pos, 1) = $n;
		$nins -= $ins;
		$tins++;
	}
	print "  Added $tins insertions, GC=$gcl\n" if(defined $debug);
	return $seq;	
}

# addTransitions => do transitions (A<=>G, T<=>C)
sub addTransitions {
	my $seq  = shift @_;
	my $nsit = shift @_;
	my $gcl  = shift @_;
	my $eval = shift @_; $eval ||= 0;
	my $tsit = 0;
	my $skip = 0;
	my @pos  = selPosition($seq, $gcl);
	while ($nsit > 0) {
		last unless(defined $pos[0]);
		my $pos = shift @pos;
		next if ($pos < $kmer);
		my $pre  = substr($seq, $pos - $kmer, $kmer + 1);
		my $post = chop $pre;
		my $old  = chop $pre;
		my $new  = '';
		
		if    ($old eq 'A') { $new = 'G'; }
		elsif ($old eq 'T') { $new = 'C'; }
		elsif ($old eq 'G') { $new = 'A'; }
		elsif ($old eq 'C') { $new = 'T'; }
		else                { next;     }
		
		my $pre_old = "$pre$old";
		my $pre_new = "$pre$new";
		my $post_old = substr("$pre_old$post", 1);
		my $post_new = substr("$pre_new$post", 1);
		
		if ($eval < $mut_cyc) {
			next unless(defined $elemk{$gcl}{$pre_old}  and defined $elemk{$gcl}{$pre_new});
			next unless(defined $elemk{$gcl}{$post_old} and defined $elemk{$gcl}{$post_new});
			next if($elemk{$gcl}{$pre_old}  >= $elemk{$gcl}{ $pre_new} + 1e-15);
			next if($elemk{$gcl}{$post_old} >= $elemk{$gcl}{$post_new} + 1e-15);
		}
		
		substr($seq, $pos - 1 , 1) = $new;
		$nsit--;
		$tsit++;
	}
	$skip = $nsit;
	print "  Added $tsit transitions ($skip skipped, GC=$gcl)\n" if(defined $debug);
	if ($skip > 0) {
		$gcl = newGC();
		$seq = addTransitions($seq, $skip, $gcl, ++$eval);
	}
	return $seq;
}

# addTransversions => do transversions (A<=>T, C<=>G)
sub addTransversions {
	my $seq  = shift @_;
	my $nver = shift @_;
	my $gcl  = shift @_;
	my $eval = shift @_; $eval ||= 0;
	my $tver = 0;
	my $skip = 0;
	my @pos  = selPosition($seq, $gcl);
	
	while ($nver > 0) {
		last unless(defined $pos[0]);
		my $pos = shift @pos;
		next if ($pos < $kmer);
		my $pre  = substr($seq, $pos - $kmer, $kmer + 1);
		my $post = chop $pre;
		my $old  = chop $pre;
		my $new  = '';
		
		if    ($old eq 'A') { $new = 'T'; }
		elsif ($old eq 'T') { $new = 'A'; }
		elsif ($old eq 'G') { $new = 'C'; }
		elsif ($old eq 'C') { $new = 'G'; }
		else                { next;       }
		
		my $pre_old = "$pre$old";
		my $pre_new = "$pre$new";
		my $post_old = substr("$pre_old$post", 1);
		my $post_new = substr("$pre_new$post", 1); 
		
		if ($eval < $mut_cyc) {
			next unless(defined $elemk{$gcl}{$pre_old}  and defined $elemk{$gcl}{$pre_new});
			next unless(defined $elemk{$gcl}{$post_old} and defined $elemk{$gcl}{$post_new});
			next if($elemk{$gcl}{$pre_old}  >= $elemk{$gcl}{ $pre_new} + 1e-15);
			next if($elemk{$gcl}{$post_old} >= $elemk{$gcl}{$post_new} + 1e-15);
		}

		substr($seq, $pos - 1, 1) = $new;
		$nver--;
		$tver++;
	}
	$skip = $nver;
	print "  Added $tver transversions ($skip skipped, GC=$gcl)\n" if(defined $debug);
	if ($skip > 0) {
		$gcl = newGC();
		$seq = addTransversions($seq, $skip, $gcl, ++$eval);
	}
	return $seq;
}

# calGC => calculate the GC content
sub calcGC {
	my $seq = shift @_;
	$seq =~ s/[^ACGTacgt]//g;
	my $tot = length $seq;
	my $ngc = $seq =~ tr/GCgc//;
	my $pgc  = int($ngc * 100 / $tot);
	my $new_gc = 40;
	
	if    ($pgc <  5) { $new_gc =   5; }
	elsif ($pgc < 10) { $new_gc =  10; }
	elsif ($pgc < 15) { $new_gc =  15; }
	elsif ($pgc < 20) { $new_gc =  20; }
	elsif ($pgc < 25) { $new_gc =  25; }
	elsif ($pgc < 30) { $new_gc =  30; }
	elsif ($pgc < 35) { $new_gc =  35; }
	elsif ($pgc < 40) { $new_gc =  40; }
	elsif ($pgc < 45) { $new_gc =  45; }
	elsif ($pgc < 50) { $new_gc =  50; }
	elsif ($pgc < 55) { $new_gc =  55; }
	elsif ($pgc < 60) { $new_gc =  60; }
	elsif ($pgc < 65) { $new_gc =  65; }
	elsif ($pgc < 70) { $new_gc =  70; }
	elsif ($pgc < 75) { $new_gc =  75; }
	elsif ($pgc < 80) { $new_gc =  80; }
	elsif ($pgc < 85) { $new_gc =  85; }
	elsif ($pgc < 90) { $new_gc =  90; }
	elsif ($pgc < 95) { $new_gc =  95; }
	else              { $new_gc = 100; }
	
	return $new_gc;
}

# newGC => obtain a new GC based on probabilities
sub newGC {
	my $gc = $classgc[0];
	if ($#classgc > 1) {
		my $dice   = rand();
		my $p      = 0;
		foreach my $class (@classgc) {
			my $q = $p + $gc{$class};
			$gc   = $class if($dice >= $p);
			last if($dice >= $p and $dice <= $q);
			$p    = $q;
		}
	}
	return $gc;
}

# transGC => get a new GC change based on probabilities
sub transGC {
    my $old_gc = shift @_;
    my $new_gc = $old_gc;
    my $dice   = rand();
    my $p      = 0;
    foreach my $gc (keys %{ $gct{$old_gc} }) {
        my $q   = $p + $gct{$old_gc}{$gc};
		$new_gc = $gc if($dice >= $p);
		last if($dice >= $p and $dice <= $q);
		$p      = $q;
	}
    return $new_gc;
}

# createSeq => first step to create a new sequence (major loop)
sub createSeq {
	my $k     = shift @_;
	my $gc    = shift @_;
	my $len   = shift @_;
	my $win   = shift @_;
	my $seq   = shift @_;
	print "creating new sequence\n" if (defined $debug);
	
	for (my $i = length $seq; $i <= $len; $i += $win) {
	    warn "    $i fragment, GC=$gc\n" if (defined $debug);
		my $seed   = substr($seq, 1 - $k);
		my $subseq = createSubSeq($k, $gc, $win - $k + 1, $seed);
		$seq      .= $subseq;
		$gc        = transGC($gc);
	}
	return $seq;
}

# createSeq => second step to create a new sequence (minor loop)
sub createSubSeq {
	my $k = shift @_;
	my $g = shift @_; 
	my $w = shift @_; 
	my $s = shift @_;

	# Extent to the window
	for (my $i = length $s; $i <= $w; $i++) {
		my $seed = substr ($s, 1 - $k);
		my $dice = rand();
		my $n    = $dna[$#dna];
		my $p    =  0;
		foreach my $b (@dna) {
		    my $q = $p + $elemk{$g}{"$seed$b"};
			$n    = $b if ($dice >= $p);
			last if($dice >= $p and $dice <= $q);
			$p    = $q;
		}
		$s .= $n;
	}
	return $s;
}

# insertElements => insert the elements
sub insertElements {
    warn "inserting elements\n" if (defined $debug);
	my $s    = shift @_;
	my %pos  = randSel((length $s) - $win, $nrep + $nsim);
	my @ins  = ();
	for (my $i = 0; $i <= $nrep; $i++) { push @ins, 'rep'; }
	for (my $i = 0; $i <= $nsim; $i++) { push @ins, 'sim'; }
	@ins = shuffle(@ins);
	foreach my $pos (keys %pos) {
	    warn "Insertion position:$pos\n" if (defined $debug);
	    my $zone = substr ($s, $pos, $win);
	    my $gc   = calcGC($zone);
	    
	    warn "  region GC=$gc\n" if (defined $debug);
	    my $ins  = shift @ins;
	    my $new  = '';
	    my $seq  = '';
	    my @bag  = ();
	    if ($ins eq 'sim') {
	        @bag = @{ $simple{$gc} };
	        $new = $bag[int(rand @bag)];
	        warn "  $new\n" if (defined $debug);
	        $seq = evolveSimple($new, $gc);
	    }
	    else {
	        @bag = @{ $repeat{$gc} };
	        $new = $bag[int(rand @bag)];
	        warn "  $new\n" if (defined $debug);
	        $seq = evolveRepeat($new, $gc);
	        next if ($seq eq 'BAD');
	    }
	    
	    $seq = lc $seq if (defined $mask);
	    
	    if ($seq =~ m/X/i) {
	        my @frag = split (/X/i, $seq);
	        my @pos  = ();
	        foreach my $frag (@frag) {
	            push @pos, $pos;
	            substr ($s, $pos, length $frag) = $frag;
	            $pos += (length $frag) + int(rand (length $frag)) + int(rand (length $frag));
	        }
	        $pos = join ":", @pos;
	    }
	    else {
	        substr ($s, $pos, length $seq) = $seq;
	    }
	    push @inserts, "$pos\t$new";
	}
	return $s;
}

# evolveSimple => return the evolved repeat
sub evolveSimple {
    my $sim   = shift @_;
    my $gc    = shift @_;
    my $seq   = '';
    my ($lab, $seed, $exp, $div, $indel) = split (/:/, $sim);
    $seq      = $seed x (int($exp) + 1);
    $seq      = substr ($seq, 0, int($exp * length $seed));
    my $mut   = int($div * (length $seq) / 100);
    my $nsit  = int($mut / 2);
    my $nver  = $mut - $nsit;
    my $nid   = int($indel * (length $seq) / 100);
    my $ndel  = int(rand $nid);
    my $nins  = $nid - $ndel;
    $seq = addDeletions(    $seq, $ndel, $gc) if($ndel > 0);
    $seq = addTransitions(  $seq, $nsit, $gc) if($nsit > 0);
    $seq = addTransversions($seq, $nver, $gc) if($nver > 0);
    $seq = addInsertions(   $seq, $nins, $gc) if($nins > 0);
    return $seq;
}

# evolveRepeat => return the evolved repeat
sub evolveRepeat {
    my $rep   = shift @_;
    my $gc    = shift @_;
    my $seq   = '';
    my ($type, $fam, $dir, $div, $ins, $del, $ini, $end) = split (/:/, $rep);
    my ($mut, $nins, $ndel, $nsit, $nver, $cseq);
    
    unless (defined $rep_seq{$type}) {
        warn "sequence for $type ($fam) not found!";
        return "BAD";
    }
    
    if ($rep  =~ m/;/) {
        my $sseq = '';
        my @frag  = split (/;/, $rep);
        my $first = shift @frag;
        ($type, $fam, $dir, $div, $ins, $del, $ini, $end) = split (/:/, $first);
        $cseq  = $rep_seq{$type};
        
        if (length $cseq < $ini) {
            warn "sequence for $type ($fam) too short!";
            $ini = 1;
        }
        
        if (length $cseq < $end) {
            $end = length $cseq;
        }
        
        $sseq = substr ($cseq, $ini - 1, $end - $ini);
        $sseq = revcomp($sseq) if ($dir eq '-');
        $mut  = int($div * (length $sseq) / 100);
        $nsit = int($mut / 2);
        $nver = $mut - $nsit;
        $nins = int($ins * (length $sseq) / 100);
        $ndel = int($ins * (length $sseq) / 100);
        $sseq = addDeletions(    $sseq, $ndel, $gc) if($ndel > 0);
        $sseq = addTransitions(  $sseq, $nsit, $gc) if($nsit > 0);
        $sseq = addTransversions($sseq, $nver, $gc) if($nver > 0);
        $sseq = addInsertions(   $sseq, $nins, $gc) if($nins > 0);
        $seq  = $sseq;
        foreach my $frag (@frag) {
            ($div, $ins, $del, $ini, $end) = split (/:/, $frag);
            if (length $cseq < $ini) {
                warn "sequence for $type ($fam) too short!";
                $ini = 1;
            }
        
            if (length $cseq < $end) {
                $end = length $cseq;
            }
            
            $sseq  = substr ($cseq, $ini - 1, $end - $ini);
            $sseq  = revcomp($sseq) if ($dir eq '-');
            $mut  = int($div * (length $sseq) / 100);
            $nsit = int($mut / 2);
            $nver = $mut - $nsit;
            $nins = int($ins * (length $sseq) / 100);
            $ndel = int($ins * (length $sseq) / 100);
            $sseq  = addDeletions(    $sseq, $ndel, $gc) if($ndel > 0);
            $sseq  = addTransitions(  $sseq, $nsit, $gc) if($nsit > 0);
            $sseq  = addTransversions($sseq, $nver, $gc) if($nver > 0);
            $sseq  = addInsertions(   $sseq, $nins, $gc) if($nins > 0);
            $seq  .= "X$sseq";
        }
    }
    else {
        $cseq = $rep_seq{$type};
        if (length $cseq < $ini) {
            warn "sequence for $type ($fam) too short!";
            $ini = 1;
        }
        
        if (length $cseq < $end) {
            $end = length $cseq;
        }
        
        $seq  = substr ($cseq, $ini - 1, $end - $ini);
        $seq  = revcomp($seq) if ($dir eq '-');
        $mut  = int($div * (length $seq) / 100);
        $nsit = int($mut / 2);
        $nver = $mut - $nsit;
        $nins = int($ins * (length $seq) / 100);
        $ndel = int($ins * (length $seq) / 100);
        $seq  = addDeletions(    $seq, $ndel, $gc) if($ndel > 0);
        $seq  = addTransitions(  $seq, $nsit, $gc) if($nsit > 0);
        $seq  = addTransversions($seq, $nver, $gc) if($nver > 0);
        $seq  = addInsertions(   $seq, $nins, $gc) if($nins > 0);
    }
    return $seq;
}

# randSel => select a random numbers in a finite range
sub randSel {
	my $total = shift @_;
	my $want  = shift @_;
	my %select = ();
	for (my $i = 0; $i <= $want; $i++) {
		my $num = int(rand $total);
		if (defined $select{$num}) {
			$i--;
		}
		else {
			$select{$num} = 1;
		}
	}
	return %select;
}

# errorExit => print an error message and finish the program (return signal: 1)
sub errorExit {
	my $err_message = shift @_;
	print "ABORTED: $err_message\n";
	exit 1;
}


