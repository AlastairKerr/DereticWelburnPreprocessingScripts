#!/bin/bash 

#The program is best used in conjunction with gnu parallel
# parallel -a listOfIDs ./run_iupred.sh {} IndexedFastaFile > output

#Config
IUPRED=./iupred   # Location of iupred program 
IUPred_PATH=$PWD   # Set the environmental var for use in the program
ParseScript=./parse_iupred.pl #Location of perl script ( Ensembl Specific)  
SAMTOOLS=samtools  # Location of samtools program
tmpfile=tmp.$$  #Set a temp file with a semi random name 


#Arguments
ID=$1      # ID to extract
FA=$2      # Indexed fastA protein sequences


$SAMTOOLS faidx $FA $ID > $tmpfile  # Extract just the protein 
$IUPRED $tmpfile long | $ParseScript  # To STDOUT 
rm $tmpfile #clean up
