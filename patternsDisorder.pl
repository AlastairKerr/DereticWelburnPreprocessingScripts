#!/usr/bin/perl -w

use Getopt::Long;

my %opt;



#Define the protein FastaA files for each Species
%opt = ( 
    lengths => "lengths"
    );


GetOptions(\%opt, 
	   "lengths=s",
	   "uiphred=s", "fuzzpro=s"
	   );

unless($opt{uiphred} and $opt{fuzzpro}){
    print "Usage: patternsDisorder.pl \\\n", 
    "\t--length lengthFile \\\n",
    "\t--uiphred parseduiphredFile \\\n",
    "\t--fuzzpro parsedfuzzproFile\n";
    exit;
}
       
open(L,$opt{lengths}) or die "cannot open lengths file $opt{lengths}";
while(<L>){
    chomp;    
    my ($id, $l) = split;
    $len{$id} = $l;
}

open(D, $opt{uiphred}) or die "cannot open uiphred file $opt{uiphred}"; # from uiphred parse 
while(<D>){
    my ($id, $pos, $val) = split;
    next if $val < 0.5;
    $dis{$id}{$pos} = $val;
}




#load the output from fuzzpro
open(F, $opt{fuzzpro}) or die "cannot open fuzzpro file $opt{fuzzpro}";
while(<F>){
    chomp;
    my ($pid, $s, $e, $pattern, undef, $seq) = split;
    $h_fuzz{$pid}||=0;  
    $h_fuzz{$pid}++;
    push(@{$h_fz_start{$pid}},$s);
}

foreach my $pid (keys %h_fz_start){
    my @starts = uniq(@{$h_fz_start{$pid}});
    
    my $i = 300; # Code changed, was going to add a loop here for various values for $i
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
	    push(@{$cloc{$i}}, $starts[$j]);
	  L:    foreach my $k ($j+1..$#starts){
	      if($starts[$k] - $starts[$j] < $i){
		  $c{$i}++;		    
		  push(@{$cloc{$i}}, $starts[$k]);
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
    #Count the dis stuff
    @st = sort {$a <=> $b } @{$mloc{$i}};
    my $a = $st[0]-50;
    my $b = $st[$#st]+50;
    my  $c = 0;
    #if near N terminus move the box to the right
    if ($a < 1){
	$b +=  1 - $a;
	$a = 1;
    }	
    #if near C terminus move box to left 
    if ($b > $len{$pid}){
	$a -= $b - $len{$pid};
	$b = $len{$pid};
    }

#count number of dis residues in box    
    foreach my $i ($a..$b){
	$c++ if defined ($dis{$pid}{$i});
    }
    
    my $pc = (100 * $c) / ($b - $a); 
    print join("\t", $pid, $i, $m{$i}, 
	       join(",", @st ), #pattern starts
	       $a, $b, $c, int($pc) # region examined for disorder
	), "\n";
    
}



sub uniq{
    my %hash;
    $hash{$_} = 1 foreach (@_);
    return (sort {$a <=> $b} keys %hash);
}
