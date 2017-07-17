
## Reference files

### Human Orthologs and GO terms 

Ensembl Biomart was used to extract a table of orthologs pairs of human genes to chicken and zebrafish ("HumanChickenZebrafish.ortho") as well as a table of human genes that mapped to our GO terms of interest ("human_go.txt"). 

Peptide sequences for each of the three species were downloaded from the Ensembl ftp site and indexed using "samtools faidx".

"lengths" : length of each portein in any dataset.
"pattern.fz" : List of the patterns used for searching.


### Mitosis Genes
Reactome as used to extracted genes involved in mitotis pathways. 


## Computation
"fuzzpro" from the EMBOSS suite was used to find matches to each pattern in each of the 3 peptide datasets.
"tabFuzzpro.pl" tabulates the fuzzpro output (expects Ensembl identifiers).
"regionAln.pl" uses the "needle" program from the EMBOSS package to calculate the best scoring local alignment between species. Uses input from "tabFuzzpro.pl". 

"run_iupred.sh" usually used with gnu parallel e.g. "parallel -a listOfIDs ./run_iupred.sh {} IndexedFastaFile > output"
"parse_iupred.pl" tabulates output from the iupred program for amino acids with values over 0.5.
"patternsDisorder.pl" counts the number of disordered peptides within the pattern regions.  Uses input from "parse_iupred.pl" and "tabFuzzpro.pl"
"mergeSpeciesRegions.pl" merges each of the 3 output files (one per species) into a single table. 

"mergeAlnRegions.R" takes the output of "regionAln.pl", and "mergeSpeciesRegions.pl" and saves an Rdata file to load into the shiny app. 

"countPatterns.pl" used to count the maximum number of patterns that are found within certain window size. Generates the ".count" files used in the shiny app.




