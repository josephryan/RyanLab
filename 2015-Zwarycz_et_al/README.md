# Scripts/ - Scripts from Zwarycz et al, 2015

note: some of these scripts require the JFR::Fasta perl module available here:
https://github.com/josephryan/JFR-PerlModules

## hd60.hmm

This is a hidden markov model generated from all homeodomains in HomeoDB (http://homeodb.zoo.ox.ac.uk/) that are 60 amino acids in length. This is 1,563 homeoboxes from Human, Mouse, Chicken, Frog, Zebrafish, Amphioxus, Fruitfly, Beetla, Honeybee, and Nematode downloaded in July 2014.  The default homeodomain HMM from PFAM is only 57 amino acids.

## branch_lengths_filter.pl

a script used to print out taxa on a tree in Newick format that occur on a terminal branch longer than 2.3.  It also prints the terminal branch lengths.

## hmmsearch_blast_combo.pl

Accepts the output of hmmsearch and the amino acid file that was used in the hmmsearch. It returns a FASTA file with all sequences except those identified in the hmmsearch.

## parse_and_reblast_w_alignments.pl

Accepts a BLASTP output (tabbed format) and the FASTA file that was used in the original hmmsearch as well as the E-Value cutoff, a director for the output, and a FASTA file used as the BLASTP database 

## remove_gaps_from_hmmsearch_results.pl

When running hmmsearch with the -A option, an alignment is output in stockholm format. After this file is converted to FASTA format, this script will remove positions in the alignment that don't correspond with the hmm (gaps and lowercase letters).

## CSV Files (CSV_files/)

# comma-separated files that include accessions names, family, class/subclass, and accession number for each gene used in the analyses

* Capitella_teleta.csv
* Crassostrea_gigas.csv
* Eisenia_fetida.csv
* Helobdella_robusta.csv
* Lottia_gigantea.csv

## Alignments (Alignments/)

# all alignments used in this study in PHYLIP format

* ALL.phy
* HOXL.phy
* LIM.phy
* NKL.phy
* OTHER.phy
* POU.phy
* PRD.phy
* SINE.phy
* TALE.phy

## Tree Files (Trees/)

# Treefiles in Newick or Nexus format for all trees generated in this study. 

* ALL.ML.tre
* HOXL.BAYES.tre
* HOXL.ML.tre
* HOXL.consense.tre
* LIM.BAYES.tre
* LIM.ML.tre
* LIM.consense.tre
* NKL.BAYES.tre
* NKL.ML.tre
* NKL.consense.tre
* OTHER.BAYES.tre
* OTHER.ML.tre
* OTHER.consense.tre
* POU.BAYES.tre
* POU.ML.tre
* POU.consense.tre
* PRD.BAYES.tre
* PRD.ML.tre
* PRD.consense.tre
* SINE.BAYES.tre
* SINE.ML.tre
* SINE.consense.tre
* TALE.BAYES.tre
* TALE.ML.tre
* TALE.consense.tre


