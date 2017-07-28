# PLANNED ANALYSES FOR TESTING CTENOPHORE PHYLOGENY AND SIGNALS OF DIVERGENCE  
 Principle Investigator: Joseph Ryan, Steven Haddock  
 Support Personnel: Melissa DeBiasse  
 Draft or Version Number: v.1.0  
 
## 1 INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_  

Ctenophore species exist across a wide depth gradient and therefore experience a range of environmental conditions relating to pressure, temperature, salinity, and oxygenation. Previous ctenophore phylogenies have been based on a small number of loci and/or incomplete taxon sampling and some nodes in the phylogeny are poorly supported.  

### 1.2 _Rationale_  

Over evolutionary time, many marine organisms have transitioned their home ranges to and from the deep sea despite the tremendous differences between these two habitats.  Such habitat shifts required dramatic genetic and physiological changes to these animal lineages over time. Comparisons of sequences in an accurate phylogenetic framework will lead to identification of some of the most important genetic changes that drove these transitions.  

### 1.3 _Objectives_  

The overall objective is to construct an accurate phylogeny of the relationships within Ctenophora and identify the genetic events that underlie physiological tolerances and adaptations to the unique challenges of the deep sea. Specifically, we aim to estimate a species phylogeny for 35 ctenophore species using transcriptome data using both concatenation- and coalescent-based species tree estimation approaches. We will generate a dataset of orthologous loci and use these data to infer species relationships and test for signals of convergence associated with depth in gene sequences. We will also perform ancestral state reconstructions to identify depth transitions and test for evidence of positive selection on these nodes.  

## 2 STUDY DESIGN AND ENDPOINTS  

#### 2.1 Translate ctenophore nucleotide transcriptome sequences into amino acid sequences with TransDecoder v3.0.0. We set the –m flag to 50 and used the results from blast and hmmscan searches to inform the final TransDecoder prediction step.  

```
TransDecoder.LongOrfs -t [transcriptome_file] -m 50  
```

   ```
blastp -query longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > outfile.blastp.out  
```

   ```
hmmscan --cpu 1 --domtblout outfile.domtblout Pfam-A.hmm longest_orfs.pep  
```

   ```
TransDecoder.Predict -t [transcriptome_file] --retain_pfam_hits outfile.domtblout --retain_blastp_hits outfile.blastp.out
```

#### 2.2 We used the program [Alien Index](https://github.com/josephryan/alien_index) to remove any contaminating, non-metazoan sequences. The blast database (ai.fa) used in the first step is available in this repository. The metazoan and non-metazoan sequences used to create the ai.fa database, and instructions for how to create a custom database, can be found [here](https://github.com/josephryan/alien_index). 

```
blastp -query [infile.pep.fa] -db ai.fa -outfmt 6 -max_target_seqs 1000 -seg yes -evalue 0.001 -out [file.out] > file.std 2> file.err
```

```
./alien_index --blast=[file_ai.out] --alien_pattern=ALIEN [out.alien_index] 2> file.std 3> file.err 


```
remove_aliens.pl [out.alien_index] [original_transcriptome.fa] > [filtered_transcriptome.fa] 2> file.std 3> file.err

```

#### 2.3 We identified orthogroups across ctenophore transcriptomes in OrthoFinder v1.1.4.  

```
orthofinder -f [dir_w_protein_fast_files] -op  
```

```
blastp -outfmt 6 -evalue 0.001 -query [renamed_fasta_file_w_all_seqs] -db BlastDB -out outfile.txt > file.std 2> file.err
```

```
orthofinder -b [dir_w_blast_results]
```

```
python trees_from_MSA.py [dir_w_orthofinder_results]
```

#### 2.4 Our data set contains two individuals from the species *Nepheloctena* ‘red’. If *Nepheloctena* ‘red’ sp. 1 and 2 were present in an orthogroup, we used the script ```condense_nephred.pl``` to remove *Nepheloctena* ‘red’ sp. 2. If *Nepheloctena* ‘red’ sp. 2 was present and *Nepheloctena* ‘red’ sp. 1 was absent, we retained *Nepheloctena* ‘red’ sp. 2. The script is available in the scripts directory in this repository.  

#### 2.5 Generate single copy orthogroups. First, the script ```filter_ogs_write_scripts.pl``` (available in this repository) retains orthogroup fasta files that contain a user-specified minimum number of taxa (for this project 28 species, 80% of the total 35) and only one sequence per species, except for *Mertensia ovum* (which was has a disproportionate number of isoforms due to an very deep sequencing). Lastly, ```filter_ogs_write_scripts.pl``` automates the following processes: 

2.5.1 sequences within each orthogroup are aligned using Mafft v7.309 

```mafft-linsi --localpair --maxiterate 1000 --thread 20 [infile]```

2.5.2 alignments are refined using Gblockswrapper v0.03 (https://goo.gl/fDjan6)

```Gblockswrapper [infile.mafft] > outfile.mafft-gb```

2.5.3 Gblockswrapper sometimes leaves blank sequences that cause issues downstream of this pipeline; the ```remove_empty_seqs``` script from this repository, removes empty sequencess and spaces from sequence lines. 

```remove_empty_seqs [outfile.mafft-gb]```

2.5.4 Maximum-likelihood orthogroup gene trees are estimated in IQTree v1.5.5 

```iqtree-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m LG -pre [output prefix]```

2.5.5 orthogroups with multiple *M. ovum* sequences are pruned in PhyloTreePruner v1.0 

```java PhyloTreePruner [infile.tree] 28 [infile.align] 0.5 u```


#### 2.6 Concatenate 944 single-copy loci filtered from step 5 to create a matrix and partition file for use in downstream phylogenomic analyses using ```fasta2phylomatrix``` (available in the scripts directory of this repository). Definition lines in each fasta file were edited (```perl -pi.orig -e 's/\|.*$//;' *.fa```) prior to running fasta2phylomatrix.  

#### 2.7 Estimate species phylogeny using concatenated and coalescent gene tree/species tree methods.  

2.7.1 Concatenated matrix, Maximum likelihood: estimate a bootstrapped (1000 ultrafast replicates) species phylogeny in IQtree v1.5.5 using the concatenated dataset. We will use the flag -m MFP+MERGE to find best partition scheme incl. FreeRate heterogeneity and estimate the tree.

```
iqtree-omp –s [infile] –pre [prefix_for_outfiles] –nt [# of cores] –q [partition file] –m MFP+MERGE –bb 1000 –bspec GENESITE
```

2.7.2 Concatenated matrix, Bayesian inference: estimate species phylogeny in PhyloBayes-MPI v1.7 using the concatenated dataset. If PhyloBayes is not close to convergence after 1 month runtime, we will use the jackknife approach described in Simion et al. 2017.  

```
mpirun -n [# cores] pb_mpi -d [infile.phy] -cat -gtr chain1 > chain1.out 2> chain1.err
mpirun -n [# cores] pb_mpi -d [infile.phy] -cat -gtr chain2 > chain2.out 2> chain2.err
bpcomp -x [burnin] [sample_every_x_number_of_trees] <chain1> <chain2>
```

2.7.3 Coalescent-based phylogeny: estimate the species phylogeny using ASTRAL-II v4.11.1 and ASTRID v1.4. 

> i) Generate individual maximum-likelihood gene trees in IQtree. 

```
iqtree-omp –s [infile] –pre [prefix_for_outfiles] –nt [# of cores] –q [partition file] –m MFP+MERGE –bb 1000 –bspec GENESITE
```

> ii) ASTRAL-II constrains the search space to those species trees that derive their bipartitions from the input gene trees

```
java -jar astral.jar -i [gene_trees_file] -o [output_file] > file.stdout 2> file.err 
```

> iii) ASTRID uses a distance matrix generated from the input gene trees to estimate the species tree and is robust to missing data

```
ASTRID –i [infile] –o [outfile] –m bionj > file.stdout 2> file.err
```

> iv) Compute branch support using local posterior probabilities.  

#### 2.8 Infer ancestral states for depth across the ctenophore phylogeny to identify depth transitions (shallow to deep and deep to shallow) and the depth state of the most recent common ctenophore ancestor

2.8.1 We will use SIMMAP to conduct character-mapping analyses under the explicit statistical models for character evolution described in SIMMAP implemented in phytools. SIMMAP uses stochastic mutational mapping to simulate the evolution of characters on a posterior distribution of trees, resulting in estimates of posterior probability (PP) for the presence or absence of each trait (i.e., depth) at each node.  

#### 2.9 Identify lineages, genes, and sites under strong positive selection using the above orthogroups and ancestral state results.

2.9.1 Convert aligned protein sequences to nucleotide sequences with PAL2NAL v14. We used the wrapper script ```pal2nal_wrapper.pl``` (available in this repository) to match sequences from the nucleotide transcriptomes to the corresponding orthogrogroup amino acid sequences and execute PAL2NAL. 

2.9.2 use aBSREL (Smith et al, 2015) from the HyPhy package v2.2.4 to rank lineages in terms of episodic diversification along each branch of the ctenophore phylogeny.
We will use batch_files/ABSREL.bf in this repository.

2.9.3 Use BUSTED (Murrell et al, 2015) from the HyPhy package v2.2.4 to identify gene-wide identification of episodic selection across our dataset.
We will use batch file in batch_files/BUSTED.bf in this repository.

2.9.4 Use FUBAR (Murrell et al, 2013) from the HyPhy package v2.2.4 to identify specific sites under positive selection.
We will use batch_files/FUBAR.bf in this repository.

2.9.5 Use RELAX (Wertheim et al, 2015) from the HyPhy package v2.2.4 to test if the strength of selection has been relaxed or intensified along a set of branches identified a priori according to the ancestral state reconstruction [section 8].
We will use batch_files/RELAX.bf in this repository.  

#### 2.10 We will test for convergence at the genic level using the SOWH test implemented in the program SOWHAT v0.36. We will create a constraint tree such that the difference between total habitat depth is maximized between two clades. If there is an even number of taxa, there will be the same number of taxa in each clade; if there is an odd number, the taxa with the middle-depth will be assigned to the clade containing the closest depth to the middle-depth. We will use the SOWH test to compare the unconstrained gene tree to this “split-depths-constrained” tree, as well as to a species-topology-constrained tree. We will generate a metric for each orthogroup, which will be the SOWH “rank of test statistic” from the “depth-constrained” test, minus the SOWH “rank of test statistic” from the species tree topology test. A high metric is consistent with high levels convergence according to depth. Here are the commands:  

```
sowhat –-constraint [depth_constraint_tree] --aln [single_gene_alignment] --name [file_name] --raxml_model [RAxML_model] > file.stdout 2> file.err
```

```
sowhat –-constraint [species_constraint_tree] --aln [single_gene_alignment] --name [file_name] --raxml_model [RAxML_model] > file.stdout 2> file.err 
```

```
grep 'rank of test statistic' depth/sowhat.results.txt sp/sowhat.results.txt | perl -ne 'm/=\s+(\d+)/; push @ts, $1; $num++; if ($num == 2) { $diff = $ts[0] - $ts[1]; print "$diff\n"; }'
```

## 3 WORK COMPLETED SO FAR WITH DATES  

July XX 2017- we have completed steps 2.1-2.5  
[PUBLISHED PHYLOTOCOL VERSION 1.0 July XX, 2017]

## 4 LITERATURE REFERENCED  

Beaulieu, J. M., & O’Meara, B. C. (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic Biology, 65(4), 583-601.

Bollback, J. P. (2006). SIMMAP: Stochastic character mapping of discrete traits on phylogenies. BMC Bioinformatics, 7.

Church, S. H., Ryan, J. F., & Dunn, C. W. (2015). Automation and evaluation of the SOWH Test with SOWHAT. Systematic Biology, 64(6), 1048-1058.

Emms, D. M., & Kelly, S. (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology, 16(1), 157. 

Gblockswrapper: http://bit.ly/2svaKcR

Goolsby, E. W., Bruggeman, J., & Ané, C. (2017). Rphylopars: fast multivariate phylogenetic comparative methods for missing data and within‐species variation. Methods in Ecology and Evolution, 8(1), 22-27.

Goolsby, E. W. (2017). Rapid maximum likelihood ancestral state reconstruction of continuous characters: A rerooting‐free algorithm. Ecology and Evolution, 7(8), 2791-2797.

Kocot, K. M., Citarella, M. R., Moroz, L. L., & Halanych, K. M. (2013). PhyloTreePruner: a phylogenetic tree-based approach for selection of orthologous sequences for phylogenomics. Evolutionary Bioinformatics Online, 9, 429.

Lartillot, N., Rodrigue, N., Stubbs, D., & Richer, J. (2013). PhyloBayes MPI: phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Systematic Biology, 62(4), 611-615.

Mirarab, S., & Warnow, T. (2015). ASTRAL-II: coalescent-based species tree estimation with many hundreds of taxa and thousands of genes. Bioinformatics, 31(12), i44-i52.

Murrell, B., Weaver, S., Smith, M. D., Wertheim, J. O., Murrell, S., Aylward, A., ... & Scheffler, K. (2015). Gene-wide identification of episodic selection. Molecular Biology and Evolution, 32(5), 1365-1371.

Nielsen, R. (2002). Mapping mutations on phylogenies. Systematic Biology 51, 729-739.

Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2014). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268-274.

Revell, L. J. (2012). phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, I, 217-223.

Smith, M. D., Wertheim, J. O., Weaver, S., Murrell, B., Scheffler, K., & Kosakovsky Pond, S. L. (2015). Less is more: an adaptive branch-site random effects model for efficient detection of episodic diversifying selection. Molecular Biology and Evolution, 32(5), 1342-1353.

Suyama, M., Torrents, D., & Bork, P. (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Research, 34(suppl_2), W609-W612.

TransDecoder: https://transdecoder.github.io/

Vachaspati, P., & Warnow, T. (2015). ASTRID: accurate species trees from internode distances. BMC Genomics, 16(10), S3.

Wertheim, J. O., Murrell, B., Smith, M. D., Kosakovsky Pond, S. L., & Scheffler, K. (2014). RELAX: detecting relaxed selection in a phylogenetic framework. Molecular Biology and Evolution, 32(3), 820-832.

Yamada, K. D., Tomii, K., & Katoh, K. (2016). Application of the MAFFT sequence alignment program to large data—reexamination of the usefulness of chained guide trees. Bioinformatics, 32(21), 3246-3251.

Yang, Ziheng. "PAML 4: phylogenetic analysis by maximum likelihood." Molecular Biology and Evolution 24.8 (2007): 1586-1591.

## APPENDIX

Version&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;Date&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Significant Revisions  
1.1  
1.2  
1.3  
1.4  



























