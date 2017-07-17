#!/usr/bin/perl -w


use Getopt::Long;


%opt = ( GO    => "human_go.txt",
	 reg1  => "hs.regions", 
	 reg2  => "gg.regions", 
	 reg3  => "dr.regions",
	 ortho => "HumanChickenZebrafish.ortho"
    );

GetOptions(\%opt, "GO=s", "reg1=s","reg2=s","reg3=s","ortho=s");



# Print Header information
print join("\t", qw(Name desc GID PID1 PID2 PID3
	       H_Count H_pos H_Count_Disorder H_pc_Disorder
	       G_Count G_pos G_Count_Disorder G_pc_Disorder
	       D_Count D_pos D_Count_Disorder D_pc_Disorder
	       GO)), "\n";


open(GO, $opt{"GO"}) or die "Cannot open GO file ";

while(<GO>){
    next if /Ensembl/;
    chomp;
    my ($gid, $pid, $go, $goname, $godef, $gocode, $godom) = split /\t/; 
    $pid2gid{$pid} = $gid;
#    $pidGo{$pid}{$go} = join("\t", $goname, $godef, $gocode, $godom);
#    push(@{$pid2go{$pid}},$go);
    push(@{$gid2go{$gid}},$go);
}

sub uniq{
    my @k = @_;
    my %k;
    foreach my $k (@k){
	next unless($k =~ /GO:/);
	$k{$k} = 1;
    }
    return(keys %k);
}


foreach my $file ($opt{"reg1"} , $opt{"reg2"} , $opt{"reg3"} ){
    open(F, $file) or die "cannot open file $file";
    while(<F>){
	chomp;
	my ($pid, undef, $count, $s, undef, undef, $dis, $pcdis) = split;
	$count{$pid} = $count;
	$start{$pid} = $s;
	$cdis{$pid} = $dis;
	$pdis{$pid} = $pcdis;
    }
}


open(I, $opt{ortho} ) or die "Cannot open ortholog table";
#### head -1 *ortho| tr '\t' '\n' | less -N
 #   1 Ensembl Gene ID
 #   2 Ensembl Transcript ID
 #   3 Ensembl Protein ID
 #   4 Chicken Ensembl Gene ID
 #   5 Canonical Protein or Transcript ID
 #   6 Chicken Ensembl Protein ID
 #   7 Homology Type
 #   8 Zebrafish Ensembl Gene ID
 #   9 Canonical Protein or Transcript ID
 #   10 Zebrafish Ensembl Protein ID
 #   11 Homology Type
 #   12 Associated Gene Name
 #   13 Description
    

while(<I>){
    chomp;
    my( $hs_g, $hs_p, 
	$gg_p, $g_type, 
	$dr_p, $d_type, 
	$name, $desc ) = (split(/\t/))[0,2,
				       5,6,
				       9,10,
				       11,12];
    
    $g_type||="";
    $d_type||="";

    next unless $g_type eq "ortholog_one2one";
    next unless $d_type eq "ortholog_one2one";
    
    next unless defined $count{$hs_p}; 
    next unless defined $count{$gg_p} ; 
    next unless defined $count{$dr_p} ; 
    

    print join("\t", 
	       $name,  $desc, 
	       $hs_g, $hs_p, $gg_p, $dr_p, 
	       $count{$hs_p},	$start{$hs_p},	$cdis{$hs_p},	$pdis{$hs_p},
	       $count{$gg_p},	$start{$gg_p},	$cdis{$gg_p},	$pdis{$gg_p},
	       $count{$dr_p},	$start{$dr_p},	$cdis{$dr_p},	$pdis{$dr_p},
	       join(",", &uniq(@{$gid2go{$hs_g}})),
	), "\n";
}

