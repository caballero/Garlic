#!/usr/bin/perl

=head1 NAME 

createFakeSequence.pl

=head1 DESCRIPTION

Perl script to create a sequence similar to the intergenic regions from
the model used. The models are described in parts:

a) The composition background with Markov models from 0 to 5 in fixed windows.
b) The repeats parsed from the alignements described by RepeatMasker and 
consensus bases acording to RepBase.
c) Transitional frecuencies for fixed windows in GC% classes (0-95).
d) Pseudogenes from the genome annotation.

To create the new sequence first we need to retrieve the elements to insert 
and calculate the length of the base DNA. Then, the elements are "bombarded"
in random positions and orientation.

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
my $win     =      200; # Window size for sequence GC transition
my $help    =    undef; # Variable to activate help
my $debug   =    undef; # Variable to activate verbose mode
my $baseseq =    undef; # Size of the sequence to generate
my %model   =       (); # Hash to store model parameters
my %inserts =       (); # Hash to store the insert sequences data
my %gct     =       (); # Hash for GC transitions GC(n-1) -> GC(n)
my %elemk   =       (); # Hash for kmers probabilities
my %gc      =       (); # Hash for GC content probabilities
my @classgc =       (); # Array for class GC
my %classgc =       (); # Hash for class GC
my $mingc   =        0; # Minimal GC content to use
my $maxgc   =       95; # Maximal GC content to use
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
				'debug|d'   	=> \$debug)
);

usage() if (defined $help);
usage() unless (defined $model and defined $size and defined $out);

# Loading model parameters
%model = readConfig($dir, $model);

# GC classes creation
for (my $i = $mingc; $i <= $maxgc; $i += 5) { 
	push @classgc, $i;
	$classgc{$i}++;
}


print "Generating a $size sequence with $model model, output in $out\n" if (defined $debug);

# Checking the size (conversion of symbols)
$size = checkSize($size);
errorExit("$size isn't a number") unless($size =~ /^\d+$/ and $size > 1);
my $portion = 1;
$portion = $model{'base_portion'} if (defined $model{'base_portion'});
$baseseq = $size * $portion;

# Loading background models
loadGCt($dir, $model{'model'}, $win);
print "GC transitions in windows size $win, loaded\n" if(defined $debug);

loadKmers($dir, $model{'model'}, $kmer, $win);
print "$kmer-mers in windows size $win, loaded\n" if(defined $debug);

# Number of simple repeats to use
unless (defined $nsim and $nsim >= 0) {
	$nsim = calcInsertNum($size, $model{'simple_count'} / $model{'genome_size'});
}
print "$nsim simple repeats to select\n" if(defined $debug);
if ($nsim > 0) {
#	while (1) {
		%{ $inserts{'simple'} } = ();
		my $sim_size = selectSimple($dir, $model{'model'}, $model{'simple_file'}, $model{'simple_count'}, $nsim);
#		if (($baseseq - $sim_size) > 0) {
#			$baseseq -= $sim_size;
#			last;
#		}
#		print "Bad selection, trying again\n" if (defined $debug);
#	}
print "Selected $nsim simple repeats ($sim_size bp)\n" if(defined $debug);
}

# Number of interspearsed repeats to use
unless (defined $nrep and $nrep >= 0) {
	$nrep = calcInsertNum($size, $model{'repeats_count'} / $model{'genome_size'});
}
print "$nrep interspersed repeats to select\n" if(defined $debug);

if ($nrep > 0) {
	# loading repeats consensus
	loadRepeatConsensus("$dir/repbase/RepeatMaskerLib.embl");
#	while (1) {
		%{ $inserts{'repeat'} } = ();
		my $rep_size = selectRepeat($dir, $model{'model'}, $model{'repeats_file'}, $model{'repeats_count'}, $nrep);
#		if (($baseseq - $rep_size) > 0) {
#			$baseseq -= $rep_size;
#			last;
#		}
#		print "Bad selection, trying again\n" if (defined $debug);
#	}
print "Selected $nrep interspersed repeats ($rep_size bp)\n" if(defined $debug);
}

# Number of pseudogenes to use
#unless (defined $npseudo and $npseudo >= 0) {
#	$nsim = calcInsertNum($size, $model{'pseudo_count'} / $model{'genome_size'});
#}
#if ($npseudo > 0) {
#	while (1) {
#		%{ $inserts{'pseudo'} } = ();
#		my $pseudo_size = selectPseudo($dir, $model{'pseudo_file'}, $npseudo);
#		if (($baseseq - $pseudo_size) > 0) {
#			$baseseq -= $pseudo_size;
#			last;
#		}
#	}
#}
#print "Selected $npseudo pseudogenes" if(defined $debug);

# Generation of base sequence
my $fgc    = newGC();
my @fseeds = keys %{ $elemk{$fgc} };
my $fseed  = $fseeds[int(rand @fseeds)];

$seq = createSeq($kmer, $fgc, $baseseq, $win, $fseed);
print "Base sequence generated ($baseseq bases)\n" if(defined $debug);

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

open  INS, ">$out.inserts" or errorExit("cannot open $out.inserts");
print INS  "TYPE\tID\tPOS\tTAG\tSEQ\n";
foreach my $type (keys %inserts) {
	foreach my $nid (keys %{ $inserts{$type} }) {
		my $p = $inserts{$type}{$nid}{'pos'};
		my $t = $inserts{$type}{$nid}{'tag'};
		my $s = $inserts{$type}{$nid}{'seq'};
		print INS "$type\t$nid\t$p\t$t\t$s\n";	
	}
}
close INS;

#################################################
##      S  U  B  R  O  U  T  I  N  E  S        ##
#################################################

sub usage {
print <<__HELP__
Usage: perl intergenic.pl [--help|-h] -o MODEL -l SIZE -n OUFILE [PARAMETERS] 
Parameters:
  -o --model     Model to use (like hg18, mm9, ... etc).
  -l --length    Size in bases [kb, Mb, Gb accepted].
  -n --name      Output files to create [*.fasta and *.log].
	
Optional or automatic parameters:
  -w --win       Window size for base generation profile.     Default =   200
  -k --kmer      Seed size to use [available: 1,2,3,4,5,6].   Default =     2
  -g --mingc     Minimal GC content to use [0,5,10,..,95].    Default =     0
  -c --maxgc     Maximal GC content to use [0,5,10,..,95].    Default =    95
  -r --repeats   Number of total repeats to insert.           Default =  Auto
  -s --simple    Number of total simple repeats to insert.    Default =  Auto
  -m --mask      Mask repeats in final sequence.              Default = False
  -d --debug     Verbose output for debug.                    Default = False
  -h --help      Print this screen.
	
__HELP__
;
exit 1;
}

sub formatFasta {
	my $sseq  = shift @_;
	my $col   = shift @_;
	$col ||= 80;
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

sub revcomp {
	my $rseq =  shift @_;
	$rseq    =~ tr/ACGTacgt/TGCAtgca/;
	return reverse $rseq;
}

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

sub calcInsertNum {
	my $size = shift @_;
	my $freq = shift @_;
	my $avg = int($size * $freq);
	$avg++;
	my $num = 0;
	my @op  = qw/plus minus none/;
	my $op  = $op[int(rand(@op))];
	my $dif = int(rand($avg) * int(rand(2)));
	if    ($op eq 'plus')    { $num = $avg + $dif; }
	elsif ($op eq 'minus')   { $num = $avg - $dif; }
	else                     { $num = $avg;        }
	$num = 0 if ($num < 0);
	return $num;
}

sub readConfig {
	my $path = shift @_;
	my $mod  = shift @_;
	my %mod  = ();
	open C, "$path/$mod/$mod.conf" or errorExit("cannot open $path/$mod/$mod.conf");
	while (<C>) {
		chomp;
		next if(/^#/);
		my($param,$value) = split (/=/, $_);
		$mod{$param} = $value;
	}
	close G;
	return %mod;
}

sub checkSize{
	my $sz = shift @_;
	$sz = lc $sz;
	if    ( $sz =~ /kb/ ) { $sz =~ s/kb//; $sz *= 1000;       }
	elsif ( $sz =~ /mb/ ) { $sz =~ s/mb//; $sz *= 1000000;    }
	elsif ( $sz =~ /gb/ ) { $sz =~ s/gb//; $sz *= 1000000000; }
	else                  { $sz =~ s/\D//g;                   }
	return $sz;
}

sub loadGCt {
	my $path   = shift @_;
	my $mod    = shift @_;
	my $win    = shift @_;
	my %gc_sum = ();
	open G, "$path/$mod/$mod.GCt.W$win" or errorExit("cannot open $path/$mod/$mod.GCt.W$win");
	while (<G>) {
		chomp;
		next if (/#/);
		my($pre,$post,$p, $num) = split (/\t/, $_);
		next unless($pre  >= $mingc and $pre  <= $maxgc);
		next unless($post >= $mingc and $post <= $maxgc);
		$gct{$pre}{$post} = $num;
		$gc_sum{$pre} += $num;
	}
	close G;
	foreach my $pre (keys %gct) {
		foreach my $post (keys %{ $gct{$pre} }) {
			$gct{$pre}{$post} /= $gc_sum{$pre};
		}
	}
}

sub loadKmers {
	my $path = shift @_;
	my $mod  = shift @_;
	my $kmer = shift @_;
	my $win  = shift @_;
	my $gc   = undef; 
	my $tot  = 0;
	open K, "$path/$mod/$mod.K$kmer.W$win" or errorExit("cannot open $path/$mod/$mod.K$kmer.W$win");
	while (<K>) {
		if (/#.*N=(\d+), GC=(\d+)-/) {
			$gc      = $2;
			$gc{$gc} = $1 + 1 if(defined $classgc{$gc});
			$tot    += $1 + 1 if(defined $classgc{$gc});
		}
		else {
			chomp;
			my ($b, $v, $f) = split (/\t/, $_);
			$elemk{$gc}{$b} = $f;
		}
	}
	close K;
	# Adjust the probability of GC
	foreach $gc (keys %classgc) { $gc{$gc} /= $tot; }
}

sub loadRepeatConsensus {
	my $file = shift @_;
	my $rep  = '';
	open REF, "$file" or errorExit("Cannot open $file");
	while (<REF>) {
		chomp;
		if (/^ID/) {
			my @line = split (/\s+/, $_);
			$rep = $line[1]; 
		}
		elsif (/^\s+/) {
			chomp;
			s/\s+//g;
			s/\d+//g;
			$rep_seq{$rep} .= $_; 
		}
		else {
			#Nothing
		}
	}
	close REF;
}

sub selectSimple {
	my $path   = shift @_;
	my $mod    = shift @_;
	my $file   = shift @_;
	my $total  = shift @_;
	my $wanted = shift @_;
	my %sel    = randSel($total, $wanted);
	my $line   = 1;
	my $cnt    = 0;
	my $size   = 0;
	
	open F, "$path/$mod/$file.W$win" or errorExit("cannot open $path/$mod/$file.W$win");
	while (<F>) {
		next unless(defined $sel{$line++});
		my $gc     = newGC();
		
		chomp;
		my @line   = split (/\t/, $_);
		my $tag = join ":", @line;
		print "Simple repeat selected: $tag (GC=$gc)\n" if(defined $debug);
		my $copy   = $line[5];
		my $match  = $line[7];
		my $indel  = $line[8];
		my $mer    = $line[15];
		my $pre_gc = $line[16];
		my $post_gc= $line[17];
		my $simple = $mer x int($copy);
		my $len    = length $simple;
		
		# Avoid small sequences
		unless ((length $repeat) > 4 * $kmer) {
			print "Too short sequence in repeat $tag - $repeat\n" if(defined $debug);
			$sel{$line + int(rand $max_cyc)}++;
			next;
		}
		
		# Null values in flanking sequences GC
		if ($pre_gc =~ /NA/) { 
			if ($post_gc =~ /NA/) {
				$pre_gc = $classgc[int(rand @classgc)];
				$post_gc = $pre_gc;
			}
			else {
				$pre_gc = $post_gc;
			}
		}
		else {
			if ($post_gc =~ /NA/) {
				$post_gc = $pre_gc;
			}
		}
		
		unless (($pre_gc >= $classgc[0] and $pre_gc <= $classgc[-1]) or ($post_gc >= $classgc[0] and $post_gc <= $classgc[-1])) {
			print "Flanking sequence for $tag isn't in range! skipping this\n" if(defined $debug);
			$sel{$line + 1}++;
			next;
		}
		
		print "Original: $simple\n" if (defined $debug);
		if ($match < 100) {
			$match = int($len * (100 - $match) / 100);
			$indel = int($len * $indel / 100);
			my $ndel = int(rand $indel);
			my $nins = $indel - $ndel;
			my $nsit = int($match / 2);
			my $nver = $match - $nsit;
			print "  Evolving: Deletions=$ndel, Insertions=$nins, Transitions=$nsit, Transversions=$nver\n" if (defined $debug);
			
			$simple = addDeletions(    $simple, $ndel, $gc) if($ndel > 0);
			$simple = addTransitions(  $simple, $nsit, $gc) if($nsit > 0);
			$simple = addTransversions($simple, $nver, $gc) if($nver > 0);
			$simple = addInsertions(   $simple, $nins, $gc) if($nins > 0);
			print "Evolved: $simple\n" if (defined $debug);
		}
		else {
			print "\tno mutations\n" if(defined $debug);
		}
		
		$inserts{'simple'}{$cnt}{'seq'} = $simple;
		$inserts{'simple'}{$cnt}{'tag'} = $tag;
		$inserts{'simple'}{$cnt}{'gcf'} = "$pre_gc-$post_gc";

		$cnt++;
		$size += length $simple;
	}
	close F;
	return $size;
}

sub selectRepeat {
	my $path   = shift @_;
	my $mod    = shift @_;
	my $file   = shift @_;
	my $total  = shift @_;
	my $wanted = shift @_;
	my %sel    = randSel($total, $wanted + 1);
	my $line   = 1;
	my $cnt    = 0;
	my $size   = 0;
		
	open F, "$path/$mod/$file.W$win" or errorExit("cannot open $path/$mod/$file.W$win");
	while (<F>) {
		$line++;
		next unless(defined $sel{$line});
		my $gc     = newGC();
		
		chomp;
		my @line   = split (/\t/, $_);
		my $pdiv   = $line[1];
		my $pdel   = $line[2];
		my $pins   = $line[3];
		my $dir    = $line[7];
		my $rep    = $line[8];
		my $fam    = $line[9];
		my $rini   = $line[10];
		my $rend   = $line[11];
		my $pre_gc = $line[13];
		my $post_gc= $line[14];
		
		# Null values in flanking sequences GC
		if ($pre_gc =~ /NA/) { 
			if ($post_gc =~ /NA/) {
				$pre_gc = $classgc[int(rand @classgc)];
				$post_gc = $pre_gc;
			}
			else {
				$pre_gc = $post_gc;
			}
		}
		else {
			if ($post_gc =~ /NA/) {
				$post_gc = $pre_gc;
			}
		}
		
		unless (($pre_gc >= $classgc[0] and $pre_gc <= $classgc[-1]) or ($post_gc >= $classgc[0] and $post_gc <= $classgc[-1])) {
			print "Flanking sequence for $rep isn't in range! skipping this\n" if(defined $debug);
			$sel{$line + 1}++;
			next;
		}
		
		# Direction switch
		#my @dir = qw/+ C/;
		#$dir = $dir[int(rand @dir)];
		
		# Swap coordinates if required
		if ($rini > $rend) {
			my $tmp = $rini;
			$rini   = $rend;
			$rend   = $tmp;
		}
		my $tag = join ":", @line;
		
		unless (defined $rep_seq{$rep}) {
			print "Cannot find a sequence for $rep! skipping this\n" if(defined $debug);
			$sel{$line + int(rand $max_cyc)}++;
			next;
		}
		
		my $consensus = checkBases($rep_seq{$rep});
		unless (length $consensus >= $rend) {
			print "End coordinate larger than consensus in $rep! skipping this\n" if(defined $debug);
			$sel{$line + int(rand $max_cyc)}++;
			next;

		}
		
		my $repeat = substr($consensus, $rini - 1, $rend - $rini);
		$repeat = revcomp($repeat) if ($dir eq 'C' or $dir eq '-');
		
		# If I cannot extract a valid sequence, use the next in line
		# I detected some bad sequences in the RepBase consensus
		unless ((length $repeat) > 4 * $kmer) {
			print "Too short sequence in repeat $tag - $repeat\n" if(defined $debug);
			$sel{$line + int(rand $max_cyc)}++;
			next;
		}
		
		print "Interspersed repeat selected: $tag (GC=$gc)\nOriginal: $repeat\n" if(defined $debug);
		
		my $len    = length $repeat;
		my $ndiv   = int($len * $pdiv / 100);
		my $ndel   = int($len * $pdel / 100);
		my $nins   = int($len * $pins / 100);
		my $nsit   = int($ndiv / 2);
		my $nver   = $ndiv - $nsit;

		print "  Evolving: Deletions=$ndel, Insertions=$nins, Transitions=$nsit, Transversions=$nver\n" if (defined $debug);

		$repeat = addDeletions(    $repeat, $ndel, $gc) if($ndel > 0);
		$repeat = addTransitions(  $repeat, $nsit, $gc) if($nsit > 0);
		$repeat = addTransversions($repeat, $nver, $gc) if($nver > 0);
		$repeat = addInsertions(   $repeat, $nins, $gc) if($nins > 0);
		
		$inserts{'repeat'}{$cnt}{'seq'} = $repeat;
		$inserts{'repeat'}{$cnt}{'tag'} = $tag;
		$inserts{'repeat'}{$cnt}{'gcf'} = "$pre_gc-$post_gc";
		$cnt++;
		$size += length $repeat;
		print "Evolved: $repeat\n" if(defined $debug);
	}
	close F;
	return $size;
}

sub selPosition {
	my $seq = shift @_;
	my $gc  = shift @_;
	my @pos = ();
	my %dat = ();
	for (my $i = 0; $i <= ((length $seq) - $kmer); $i++) {
		my $seed = substr($seq, $i, $kmer);
		$dat{$i} = $elemk{$gc}{$seed};
	}
	foreach my $pos (sort { $dat{$a} <=> $dat{$b} } keys %dat) {
	 	push (@pos, $pos);
	}
	return @pos;
}

sub addDeletions {
	my $seq  = shift @_;
	my $ndel = shift @_;
	my $gcl  = shift @_;
	my $eval = shift @_; $eval ||= 0;
	my $tdel = 0;
	my $skip = 0;
	my @pos  = selPosition($seq, $gcl);
	#print "PreDeletions: $seq\n" if(defined $debug);
	#print "\t\tDeleting $ndel bases in: " if(defined $debug);
	while ($ndel > 0) {
		my $bite = 1;
		last unless(defined $pos[0]);
		my $pos  = shift @pos;
		next if ($pos < $kmer);
		next if ($pos >= length $seq);
		my $pre  = substr($seq, $pos, $kmer - 1);
		my $old  = substr($seq, $pos, $kmer);
		my $new  = substr($seq, $pos + $kmer, 1);
		#print "      changing $pos: $old -> $pre$new\n" if(defined $debug);
		if ($eval < $mut_cyc) {
			next unless(defined $elemk{$gcl}{"$pre$new"} and defined $elemk{$gcl}{"$old"});
			next if($elemk{$gcl}{"$old"} >= $elemk{$gcl}{"$pre$new"} + 1e-15);
		}
		
		substr($seq, $pos, $bite) = '';
		$ndel -= $bite;
		$tdel++;
		#print " $pos" if(defined $debug);
	}
	$skip = $ndel;
	print "  Added $tdel deletions ($skip skipped, GC=$gcl)\n" if(defined $debug);
	if ($skip > 0) {
		$gcl = newGC();
		$seq = addDeletions($seq, $skip, $gcl, ++$eval);
	}
	return $seq;
}

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
		#print "      changing $pos: $seed + $n + $post\n" if(defined $debug);
		#print "\t\tAdding $ins ($n) bases in $pos position\n" if(defined $debug);
	}
	print "  Added $tins insertions, GC=$gcl\n" if(defined $debug);
	return $seq;	
}

sub addTransitions {
	my $seq  = shift @_;
	my $nsit = shift @_;
	my $gcl  = shift @_;
	my $eval = shift @_; $eval ||= 0;
	my $tsit = 0;
	my $skip = 0;
	my @pos  = selPosition($seq, $gcl);
	#print "\t\tTransitions in: " if(defined $debug);
	#print "PreTransitions: $seq\n" if(defined $debug);
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
		#print "$      changing $pos: $pre_old->$pre_new | $post_old -> $post_new\n" if (defined $debug);
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

sub addTransversions {
	my $seq  = shift @_;
	my $nver = shift @_;
	my $gcl  = shift @_;
	my $eval = shift @_; $eval ||= 0;
	my $tver = 0;
	my $skip = 0;
	my @pos  = selPosition($seq, $gcl);
	#print "\t\tTransversions in: " if(defined $debug);
	#print "PreTransversions: $seq\n" if(defined $debug);

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
		#print "$      changing $pos: $pre_old->$pre_new | $post_old -> $post_new\n" if (defined $debug);
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

sub calcGC {
	my $seq    = shift @_;
	my $tot    = length $seq;
	my $gc     = $seq =~ tr/GCgc//;
	my $new_gc = int($gc * 100 / $tot);
	#for (my $i = 0; $i <= $#classgc; $i++) {
	#	if ($classgc[$i] >= $new_gc) {
	#		$gc = $classgc[$i];
	#		last;
	#	}
	#}
	return $new_gc;
}

sub newGC {
	my $gc = $classgc[0];
	if ($#classgc > 1) {
		$gc        = $classgc[-1];
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

sub createSeq {
	my $k     = shift @_;
	my $gc    = shift @_;
	my $len   = shift @_;
	my $win   = shift @_;
	my $seq   = shift @_;
	my $tries = 0;
	my $max   = $max_cyc;
	
	for (my $i = length $seq; $i <= $len; $i += $win) {
		print "creating new sequence " if (defined $debug);
		my $seed   = substr($seq, 1 - $k);
		my $subseq = createSubSeq($k, $gc, $win - $k + 1, $seed);
		$seq      .= $subseq;
		print "param=$k, $gc, $win, $seed, GC=", calcGC($subseq), "\n" if(defined $debug);
		
		# Transition of GC
		$gc = newGC();
		#if ($#classgc > 1) {
		# 	$gc = $classgc[int(rand(@classgc))];
		#	last if (rand() <= ($gct{$old_gc}{$gc})*($gc{$gc}));
		#	if ($tries > $max) {
		#		print "Too much tries in GC transition ... using $gc\n";
		#		$tries = 0;
		#		last;
		#	}
		#	$tries++;
		#}
	}
	return $seq;
}

sub createSubSeq {
	my $k = shift @_;
	my $g = shift @_; 
	my $w = shift @_; 
	my $s = shift @_;
	my $max   = $max_cyc;

	# Extent to the window
	for (my $i = length $s; $i <= $w; $i++) {
		my $seed = substr ($s, 1 - $k);
		my $dice = rand();
		my $n    = $dna[$#dna];
		my $p    =  0;
		unless (defined $elemk{$g}{"$seed$n"}) {
			print "Bad seed ($seed) in seq ($seq)\n" if (defined $debug);
			$n = $dna[int(rand @dna)];
		}
		else {
			foreach my $b (@dna) {
				my $q = $p + $elemk{$g}{"$seed$b"};
				$n    = $b if ($dice >= $p);
				last if($dice >= $p and $dice <= $q);
				$p    = $q;
			}
		}
		$s .= $n;
	}
	return $s;
}

sub insertElements{
	my $s = shift @_;
	foreach my $type (keys %inserts) {
		foreach my $elem (keys %{ $inserts{$type} }) {
			my $insert       = $inserts{$type}{$elem}{'seq'};
			my $tag          = $inserts{$type}{$elem}{'tag'};
			$insert          = lc $insert if(defined $mask);
			my $gc_range     = $inserts{$type}{$elem}{'gcf'};
			my ($pre, $post) = split (/-/, $gc_range);
			if ($pre > $post) {
				my $tmp = $pre;
				$pre    = $post;
				$post   = $tmp;
			}
			
			my $elem_ready   = 0;
			my $tries        = 0;
			my $pos          = 0;
			my $frag_gc      = 0;
			while ($elem_ready == 0) {
				$pos      = int(rand((length $s) - $win));
				next if ($pos < ($win / 2));
				my $frag     = substr($s, $pos - ($win / 2), $win);
				$frag_gc  = calcGC($frag);
				$elem_ready = 1 if ($frag_gc >= $pre and $frag_gc <= $post);
				$inserts{$type}{$elem}{'pos'} = $pos;
				if ($tries > $max_cyc) {
					print "Too much tries in insertion, random choice\n" if (defined $debug);
					last;
				}
				$tries++;
			}
			print "Inserted $tag in $pos (GC=$frag_gc:$pre-$post)\n" if (defined $debug);
			substr($s, $pos, 1) = $insert;
		}
	}
	return $s;
}

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
			$select{$num}++;
		}
	}
	return %select;
}

sub errorExit {
	my $mess = shift @_;
	print "ABORTED: $mess\n";
	exit 1;
}


