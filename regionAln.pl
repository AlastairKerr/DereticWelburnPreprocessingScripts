#!/usr/bin/perl -w
use strict;


use Getopt::Long;

my %opt;



#Define the protein FastaA files for each Species
%opt = ( DR_fa => "Danio_rerio.GRCz10.pep.all.fa", 
	 GG_fa => "Gallus_gallus.Galgal4.pep.all.fa",
	 HS_fa => "Homo_sapiens.GRCh38.pep.all.fa",
#Extend the features of interest by this amount each side.
	 flank => 10,

#File names for the tabulated fuzzpro data
	 Hfuzz => "hs.fuzzpro",
	 Gfuzz => "gg.fuzzpro",
	 Dfuzz => "dr.fuzzpro",
	 
#Ortholog Table extracted from ensembl
	 Orthos => "HumanChickenZebrafish.ortho"
	 );    


GetOptions(\%opt, "DR_fa=s", "GG_fa=s", "HS_fa=s",
	   "flank=i",
	   "Hfuzz=s", "Gfuzz=s", "Dfuzz=s", 
	   "Orthos=s");

my @fuzzproFiles = ($opt{Hfuzz}, $opt{Gfuzz}, $opt{Dfuzz});

sub uniq{ # removes duplicates from a list 
    my @k = @_;
    my %k;
    foreach my $k (@k){
	$k{$k} = 1;
    }
    return(keys %k);
}

###
# Process parsed  fuzzpro data 

my %fz_start;
my %fz_pat;
my %fz_end;
my %fz_seq;
my %fuzz;


foreach my $fzFile (@fuzzproFiles){
    open(F, $fzFile) or die  "Cannot open Human fuzzpro file";
    

    while(<F>){
	chomp;
	my ($pid, $s, $e, $pattern, undef, $seq) = split;
	$fuzz{$pid}||=0;  
	$fuzz{$pid}++;
	push(@{$fz_start{$pid}},$s);
	push(@{$fz_end{$pid}},$e);
	push(@{$fz_pat{$pid}},$pattern);
	push(@{$fz_seq{$pid}},$seq);
    }
}

####
my %maxStart;
my %maxGg;
my %maxDr;
###
# Process ortholog data from Ensembl Biomart
open(I, $opt{Orthos}) or die "Cannot open ortholog file";

while(<I>){
    chomp;
    my( $hs_g, $hs_p, 
	$gg_p, 	$g_type, 
	$dr_p, 	$d_type, 
	$name, $desc ) 	= (split(/\t/))
	    [0,2,
	     5,6,
	     9,10,
	     11,12];
    
    $g_type||="";
    $d_type||="";


# Just select one to one orthologs to human
    next unless $g_type eq "ortholog_one2one";
    next unless $d_type eq "ortholog_one2one";

# Only select if all have a fuzzpro pattern     
    next unless defined  $fuzz{$hs_p}; 
    next unless defined  $fuzz{$gg_p} ; 
    next unless defined  $fuzz{$dr_p} ; 
    
# dereference pattern start positionsfor clarity 
    my @hstart = @{$fz_start{$hs_p}};


# cycle through each start position     
    foreach  my $h_start ( @hstart ){
	my $g_max =0;
	my $d_max =0;

# calculate a region with flanks for use in samtools  
	my $lhs = ($h_start - $opt{flank});
	$lhs = 1 if $lhs < 1;
	my $hx = $hs_p . ":" 
	    . $lhs
	    . "-" 
	    . ($h_start + 5 + $opt{flank});

#Extract  the protein sequence fragment from the master sequence 	
	my $Hstmp = "hs." . rand(1);
	system("samtools faidx $opt{HS_fa} $hx > $Hstmp");
	foreach  my $g_start (@{$fz_start{$gg_p}}){
	    my $lhs = ($g_start - $opt{flank});
	    $lhs = 1 if $lhs < 1;
	    my  $gx = $gg_p . ":" 
		. $lhs
		. "-" 
		. ($g_start + 5 + $opt{flank});
# Same for Chicken 	    
	    my $Ggtmp = "gg." . rand(1);
	    system("samtools faidx $opt{GG_fa} $gx > $Ggtmp");

	    my $tempfile = "aln." . rand(1);
	    
# Align using needle from EMBOSS 	    
	    my $temp =  system("needle asequence=$Hstmp bsequence=$Ggtmp gapopen=3 gapextend=1 outfile=$tempfile -sprotein 2> /dev/null");
	    `rm $Ggtmp`;
	    open(A, $tempfile) or die "no align for G $gx $g_start";

#Extract percentage identity 
	   # my $g_max =0;
	#    my $g_max_start; # for debugging
	    while(<A>){
		my $pc;
		next until /Ident/;
		if(/\((.*)\%/){
		    $pc = $1;
		}else{$pc=0}
#Record the max percentage identity for all pairwise combination of region in the 2 sequences 		

		if($pc > $g_max){
		    $g_max = $pc;
#		    $g_max_start = $g_start;
		    last;		
		}

	    } `rm $tempfile`;# finish reading file 
	}
	

# Now repeat the whole process for Zebrafish
# Should have use more subroutines here I know, but it works.. 
	
	
	foreach my $d_start (@{$fz_start{$dr_p}}){
	    my $lhs = ($d_start - $opt{flank});
	    $lhs = 1 if $lhs < 1 ;		       
	    my  $dx = $dr_p . ":" 
		. $lhs
		. "-" 
		. ($d_start + 5 + $opt{flank});
	    my $Drtmp = "Dr." . rand(1);
	    system("samtools faidx $opt{DR_fa} $dx > $Drtmp");
	    my $tempfile = "aln." . rand(1); # unique name for tempfile
	    my $temp = system("needle asequence=$Hstmp bsequence=$Drtmp gapopen=3 gapextend=1 outfile=$tempfile -sprotein 2> /dev/null");
	    `rm $Drtmp`; # clean up
	    open(A, $tempfile) or die "no align for D $dx $d_start";
	    while(<A>){
		my $pc;
		next until /Ident/;
		if(/\((.*)\%/){
		    $pc = $1;
		}else{$pc=0}
		if($pc > $d_max){
		    $d_max = $pc;
		    last;		
		}
	    }`rm $tempfile`; # clean up temp file 

	}
	
	$maxStart{$hs_p} ||= 0;
	$maxGg{$hs_p} ||=0;
	$maxDr{$hs_p} ||=0;
	`rm $Hstmp`; # clean up Hs fa file	

# Filter out if regions have less than 30% identity 
	next if $g_max < 30;
	next if $d_max < 30;
	#look at values pairwise 
	if(($maxGg{$hs_p} <= $g_max) and ($maxDr{$hs_p} <= $d_max)){
	    $maxStart{$hs_p} = $h_start;
	    $maxGg{$hs_p} = $g_max;
	    $maxDr{$hs_p} = $d_max;

	}
    }

# Ignore all proteins with less than 3 matching patterns 
#    next unless $h_fuzz{$hs_p} > 2; 
#    next unless $g_fuzz{$gg_p} > 2; 
#    next unless $d_fuzz{$dr_p} > 2; 
    
# Now output the data              
    print join("\t", 
	       $name,  $desc, 
	       $hs_g, $hs_p, $gg_p, $dr_p, 
	       $fuzz{$hs_p},
	       $fuzz{$gg_p},
	       $fuzz{$dr_p},	
	       $maxStart{$hs_p},
	       $maxGg{$hs_p},
	       $maxDr{$hs_p}
	), "\n";
}
