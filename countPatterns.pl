#!/bin/perl -w

open(P, "pattern.fz") or die "need pattern.fz in pwd";
while(<P>){
    chomp;
    my ($val, $pat) = split; 
    $pat{$pat} = $val;
}


#load the output from fuzzpro
open(F, shift) or die $!;
while(<F>){
    chomp;
    my ($pid, $s, $e, $pattern, undef, $seq) = split;
    $h_fuzz{$pid}||=0;  
    $h_fuzz{$pid}++;
    push(@{$h_fz_start{$pid}},$s);
#    push(@{$h_fz_end{$pid}},$e);
    $h_fz_pat{$pid}{$s} = $pat{$pattern};
 #   push(@{$h_fz_seq{$pid}},$seq);
}

foreach my $pid (keys %h_fz_start){
    my @starts = uniq(@{$h_fz_start{$pid}});
    
    foreach my $i (20, 50, 100, 300){
	my %c;
	my %m;
	my %cloc;
	my %mloc;
	$c{$i} = 1;
	$m{$i} = 1;
	#cycle through each position 
	# $starts[$j],...
	foreach my $j (0..$#starts){
	    last if $j+1 > $#starts;
	    push(@{$cloc{$i}}, $h_fz_pat{$pid}{$starts[$j]});
	L:    foreach my $k ($j+1..$#starts){
		if($starts[$k] - $starts[$j] < $i){
		    $c{$i}++;		    
		    push(@{$cloc{$i}}, $h_fz_pat{$pid}{$starts[$k]});
		}else{
		    last L;
		}
	
	    }
	    $m{$i} ||=0; # autoviv new regions
	    
	    if ($c{$i} > $m{$i}){ # store maxima
		$m{$i} =  $c{$i};
		$mloc{$i} = $cloc{$i};
	    }
	    $c{$i} =1; # reset for next loop
	    $cloc{$i} = (); # clear 
	}
	next if $m{$i} == 1;
	print join("\t", $pid, $i, $m{$i}, 
		   join(",", @{$mloc{$i}} )
#		   join(",", sort {$a <=> $b} @{$mloc{$i}} )
	    ), "\n";
    }    
}



sub uniq{
    my %hash;
    $hash{$_} = 1 foreach (@_);
    return (sort {$a <=> $b} keys %hash);
}
