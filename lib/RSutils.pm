package RSutils;
use strict;
# All methods use zero-based start, one-based end.

sub new {
	my $package = shift;
	
	my $obj = {};
	bless $obj, $package;
	return $obj;
}


sub RSload {
	#loads a range set from a file
	#the file can be gzipped
	#optional second parameter indicates the delimiter to use
	#optional third to fifth parameters indicate which columns to use for start, end and label
	my($self, $file, $delimiter, $startcol, $endcol, $labelcol) = @_;
	my(@rs, @stuff);
	$delimiter ||= "\t";
	$startcol ||= 0;
	$endcol ||= 1;
	$labelcol ||= 2;
	if ($file =~ /\.gz$/) {
		open F, "gunzip -c $file |";
	} else {
		open F, $file;
	}
	while (<F>) {
		chomp;
		@stuff = split /$delimiter/;
		push @rs, join("\t", $stuff[$startcol], $stuff[$endcol], $stuff[$labelcol]);
	}
	
	close F;
	return \@rs;
}

sub RSdump {
	#dumps the content of a range set to standard output
	#optional second parameter indicates which delimiter to use
	my($self, $RSref, $delimiter) = @_;
	$delimiter ||= "\t";
	foreach my $item (@{$RSref}) {
		if ($delimiter eq "\t") {
			print "$item\n";
		} else {
			my($start, $end, @labels) = split /\t/, $item;
			print join($delimiter, $start, $end, @labels), "\n";
		}
	}
}




sub RSsort {
	#sorts a range set by starting position and then by ending position
	#overlapping ranges are not collapsed, labels are not modified
	my($self, $RSref) = @_;
	my(@sorted, $start, $end, @starts, @ends, $i);
	
	foreach my $item (@{$RSref}) {
		($start, $end) = split /\t/, $item;
		push @starts, $start;
		push @ends, $end;
	}
	foreach $i (sort {$starts[$a]<=>$starts[$b] || $ends[$a]<=>$ends[$b]} (0..$#{$RSref})) {
	 	push @sorted, $$RSref[$i];
	}
	
	return \@sorted;
}

sub RSunion {
	#creates the union of range sets
	#can be used on just one range set to collapse overlapping ranges
	#the resulting set is returned sorted
	#repeated labels are condensed; labels are sorted per range
	my($self, @RSreflist) = @_;
	my($RSref, $item, $start, $end, @rest, %info, $laststart, $lastend, @sorted, @info);
	
	foreach $RSref (@RSreflist) {
		foreach $item (@{$RSref}) {
			($start, $end, @rest) = split /\t/, $item;
			push @{$info{$start}{$end}}, @rest;
		}
	}
	
	$laststart = "x";
	foreach $start (sort {$a<=>$b} keys %info) {
		foreach $end (sort {$a<=>$b} keys %{$info{$start}}) {
			if ($laststart eq "x") {
				$laststart = $start;
				$lastend = $end;
				@info = @{$info{$start}{$end}};
			} elsif ($start>$lastend) {
				my %labels;
				foreach (@info) {
					$labels{$_}++;
				}
				push @sorted, join("\t", $laststart, $lastend, sort keys %labels);
				$laststart = $start;
				$lastend = $end;
				@info = @{$info{$start}{$end}};
			} else {
				$lastend = $end if $end>$lastend;
				push @info, @{$info{$start}{$end}};			
			}
		}
	}
	push @sorted, join("\t", $laststart, $lastend, @info) if $laststart ne "x";
	return \@sorted;
}

sub RSintersection {
	#calculates the intersection of any number of range sets
	#if given just one set, behaves like RSunion
	#currently, it returns just the ranges, losing all extra information
	my($self, $joint, @RSreflist) = @_;
	my($RSref, @newjoint, $rangen, $lastrangen, $item, $start, $end, $rstart, $rend, $sorted);
	
	$joint = $self->RSunion($joint);
	foreach $RSref (@RSreflist) {
		@newjoint = ();
		$sorted = $self->RSunion($RSref);
		$lastrangen = 0;
		foreach $item (@{$sorted}) {
			($start, $end) = split /\t/, $item;
			for ($rangen=$lastrangen;$rangen<=$#{$joint};$rangen++) {
				($rstart, $rend) = split /\t/, $$joint[$rangen];
				last if $rstart>=$end;
				$lastrangen = $rangen;
				next if $start>$rend;
				push @newjoint, join("\t", $start<$rstart ? $rstart : $start, $end<$rend ? $end : $rend);
			}
		}
		$joint = \@newjoint;
	}
	return $joint;
}
sub RSintersection2 {
	#calculates the intersection of any number of range sets
	#assumes all range sets have already been passed by RSunion
	#currently, it returns just the ranges, losing all extra information
	my($self, $joint, @RSreflist) = @_;
	my($RSref, @newjoint, $rangen, $lastrangen, $item, $start, $end, $rstart, $rend, $sorted);
	
	#$joint = $self->RSunion($joint);
	foreach $sorted (@RSreflist) {
		@newjoint = ();
		$lastrangen = 0;
		foreach $item (@{$sorted}) {
			($start, $end) = split /\t/, $item;
			for ($rangen=$lastrangen;$rangen<=$#{$joint};$rangen++) {
				($rstart, $rend) = split /\t/, $$joint[$rangen];
				last if $rstart>=$end;
				$lastrangen = $rangen;
				next if $start>$rend;
				push @newjoint, join("\t", $start<$rstart ? $rstart : $start, $end<$rend ? $end : $rend);
			}
		}
		$joint = \@newjoint;
	}
	return $joint;
}

sub RSsubtraction {
	#calculates what's left of a range set when subtracting from it any number of range sets
	#if given just one set, behaves like RSunion
	##modified ranges retain labels - some labels may have become irrelevant
	my($self, $joint, @RSreflist) = @_;
	my($RSref, $rangen, $lastrangen, $item, $start, $end, $rstart, $rend, @replacement, $sorted, @labels);
	
	$joint = $self->RSunion($joint);
	foreach $RSref (@RSreflist) {
		$sorted = $self->RSunion($RSref);
		$lastrangen = 0;
		foreach $item (@{$sorted}) {
			($start, $end) = split /\t/, $item;
			for ($rangen=$lastrangen;$rangen<=$#{$joint};$rangen++) {
				($rstart, $rend, @labels) = split /\t/, $$joint[$rangen];
				last if $rstart>=$end;
				$lastrangen = $rangen;
				next if $start>=$rend;
				@replacement = ();
				if ($rstart<$start) {
					push @replacement, join("\t", $rstart, $rend<$start ? $rend : $start, @labels);
				}
				if ($rend>$end) {
					push @replacement, join("\t", $end, $rend, @labels);
				}
				splice @$joint, $rangen, 1, @replacement;
				$rangen--;
			}
		}
	}
	return $joint;
}

sub RSsuppression {
	#calculates what's left of a range set when removing from it any range
	#that overlaps ranges in any number of range sets
	#if given just one set, behaves like RSunion
	my($self, $joint, @RSreflist) = @_;
	my($RSref, $rangen, $lastrangen, $item, $start, $end, $rstart, $rend, $sorted);
	
	$joint = $self->RSunion($joint);
	foreach $RSref (@RSreflist) {
		$sorted = $self->RSunion($RSref);
		$lastrangen = 0;
		foreach $item (@{$sorted}) {
			($start, $end) = split /\t/, $item;
			for ($rangen=$lastrangen;$rangen<=$#{$joint};$rangen++) {
				($rstart, $rend) = split /\t/, $$joint[$rangen];
				last if $rstart>=$end;
				$lastrangen = $rangen;
				next if $start>=$rend;
				splice @$joint, $rangen, 1;
			}
		}
	}
	return $joint;
}

sub SRSsort {
	#sorts a spliced range set by starting position of the first segment and then by ending position of first segment
	#overlapping ranges are not collapsed, labels are not modified
	my($self, $RSref) = @_;
	my(@sorted, $start, $end, @starts, @ends, $i);
	
	foreach my $item (@{$RSref}) {
		($start, $end) = split /\t/, $item;
		($start) = split /,/, $start;
		($end) = split /,/, $end;
		push @starts, $start;
		push @ends, $end;
	}
	foreach $i (sort {$starts[$a]<=>$starts[$b] || $ends[$a]<=>$ends[$b]} (0..$#{$RSref})) {
	 	push @sorted, $$RSref[$i];
	}
	
	return \@sorted;
}

sub SRScluster {
	#clusters spliced ranges in one or more sets, by range overlap (e.g. cluster transcripts by exon overlap)
	my($self, @RSreflist) = @_;
	my(%info);
	foreach my $RSref (@RSreflist) {
		foreach my $item (@{$RSref}) {
			my ($starts, $ends, $id) = split /\t/, $item, 3;
			my @starts = split /,/, $starts;
			my @ends = split /,/, $ends;
			my $start = $starts[0];
			my $end = $ends[$#ends];
			
			if (defined $info{$starts}{$ends}) {
				#we already saw such an exact spliced range definition, simply append id info to the old one
				$info{$starts}{$ends}{'id'} .= ",$id";
				next;
			}
			
			#register the range
			$info{$starts}{$ends}{'start'} = $start;
			$info{$starts}{$ends}{'end'} = $end;
			$info{$starts}{$ends}{'id'} = $id;

			#compare to previous ranges...
			foreach my $pstarts (keys %info) {
				foreach my $pends (keys %{$info{$pstarts}}) {
					next if $info{$pstarts}{$pends}{'id'} eq $id;
					next if $start+1>$info{$pstarts}{$pends}{'end'};
					next if $end<$info{$pstarts}{$pends}{'start'}+1;
					#there is genomic overlap, test exons...
					my @pstarts = split /,/, $pstarts;
					my @pends = split /,/, $pends;
					foreach my $i (0..$#starts) {
						foreach my $j (0..$#pstarts) {
							last if $pstarts[$j]+1>$ends[$i];
							next if $pends[$j]<$starts[$i]+1;
							#there is exon overlap!
							my $me = join("\t", $starts, $ends);
							my $slavestarts = $pstarts;
							my $slaveends = $pends;
							my $previousMaster = $info{$slavestarts}{$slaveends}{'master'};
							last if $previousMaster eq $me;
							$info{$slavestarts}{$slaveends}{'master'} = $me;
							while ($previousMaster && ($previousMaster ne $me)) {
								($slavestarts, $slaveends) = split /\t/, $previousMaster;
								$previousMaster = $info{$slavestarts}{$slaveends}{'master'};
								$info{$slavestarts}{$slaveends}{'master'} = $me;
							}
						}
					}
				}
			}
		}
	}
	
	my %clusters;
	foreach my $starts (keys %info) {
		foreach my $ends (keys %{$info{$starts}}) {
			#identify ultimate master
			my $master = $info{$starts}{$ends}{'master'};
			my $umaster = join("\t", $starts, $ends);
			my @remaster = ($umaster);
			while ($master) {
				push @remaster, $master;
				my($mstarts, $mends) = split /\t/, $master;
				if ($master = $info{$mstarts}{$mends}{'master'}) {
					last if $master eq $umaster;
					$umaster = $master;
				}
			}
			#remaster to ultimate master and register in cluster
			foreach my $id (@remaster) {
				my($idstarts, $idends) = split /\t/, $id;
				$info{$idstarts}{$idends}{'master'} = $umaster;
				$clusters{$umaster}{$id} = 1;
			}
		}
	}
	
	#merge exons
	my @result;
	foreach my $master (keys %clusters) {
		my(@ranges, @ids);
		foreach my $id (keys %{$clusters{$master}}) {
			my($starts, $ends) = split /\t/, $id;
			my @starts = split /,/, $starts;
			my @ends = split /,/, $ends;
			foreach my $i (0..$#starts) {
				push @ranges, join("\t", $starts[$i], $ends[$i]);
			}
			push @ids, $info{$starts}{$ends}{'id'};
		}
		my $joint = $self->RSunion(\@ranges);
		my(@starts, @ends);
		foreach my $range (@{$joint}) {
			my($start, $end) = split /\t/, $range;
			push @starts, $start;
			push @ends, $end;
		}
		push @result, join("\t", join(",", @starts), join(",", @ends), join(",", @ids));
	}
	return \@result;
}


1;
