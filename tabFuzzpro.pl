#!/usr/bin/perl -w

#
# fuzzpro file as the 1st argument and outputs a tab delimeted file to STDOUT
#
#
while(<>){
    if(/Sequence:\s+([\w\d]+)/){
	$pid = $1
    }
    if(/pattern:/){
	print "$pid\t$_";
    }
}

