#!/usr/bin/perl -w 

while(<>){
    s/^ *//;
    chomp;
    my ($pos, $aa, $val) = split;
    if(/(ENS\w*P\d+)/){
	$id = $1;
    }
    /\#/  and next;
    if($val >= 0.5){
	print join("\t",  $id, $pos, $val), "\n";
    }
}
