# a place for the ryanlab to test git 

# Sea Cucumbers! 
 Principle Investigator: Joseph Ryan, Gustav Paulay 
 Support Personnel: Jessica Whelpley  
 Draft or Version Number: v.1.0  
 Date: 5 October 2017  
 Note: this document will be updated (updates will be tracked through github)
 
## 1 INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_

Sea cucumbers are the most apomorphic echinoderms: bilaterally symmetrical worms with a reduced skeleton. They have evolved bizarre specializations such as anal suspension feeding, evisceration, sticky Cuvierian tubules that entangle attackers, and a “melting” body wall. They are abundant, ubiquitous in the benthos, from poles to equator, intertidal to the deepest trenches, and include >1700 species in 25 families. They are among the most conspicuous mobile invertebrates on reefs and the deep sea and constitute the largest invertebrate fishery on Pacific islands, with stocks fully depleted throughout the tropics. A recent six gene molecular phylogeny showed substantial conflict with previous morphological-based relationships (Miller et al. 2017). Here, we present a the pipeline of our phylogenomic analysis of Holothuroidea using 9 unpublished and 16 published transcriptomes.

### 1.2 _Rationale_ 

An in-depth understanding of phylogenetic relationships is essential to understand how holothurians have adapted to the varied environments that they inhabit. Additionally, comprehension of their evolutionary relationships can directly influence fishing and conservation initiatives. The combination of a backbone phylogeny built with hundreds of genes and high-quality baits for target enrichment, will help bring phylogenetic resolution to this fascinating group of animals and provide an important set of resources for systematists to conduct low-cost phylogenetic and population sampling. 

### 1.3 _Objectives_  

The overall objective is to test the recently proposed topology (Miller et al. 2017) of higher-level relationships of holothurians. With the sea cucumber orthogroups identified through our phylogenetic analysis, and our recently published draft *Australostichopus mollis* genome (Long et al. 2016) we will design a set of target-enrichment baits for Holothuroidea. 

## 2 STUDY DESIGN AND ENDPOINTS 

#### 2.1 Search for publically available transcriptomes from SRA database (https://www.ncbi.nlm.nih.gov/sra) using the following parameters: 
```
taxa name [ORGN] AND Illumina BUTNOT genotyping BUTNOT metagenome BUTNOT mRNA
```

2.1.1 Download and split the seventeen SRR sequences using SRA toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) then rename and compress the FASTQ files. 

```
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR[number] > fqd.SRR[number].out 2> fqd.SRR[number].err &
```

```
mv SRR[number]_1.fastq SRR[number]_pass_1.fastq
```

```
pigz -9 SRR[number]_pass_1.fastq
```

#### 2.2 Use BL-FILTER on the publically available transcriptomes and our transcriptomes, part of the Agalma pipeline (Dunn et al. 2013) to trim the adapters added during Illumina RNA-Seq.

2.2.1 BL-FILTER and organizational steps for downloaded transcriptomes. Similar steps are applied to our transcriptomes, with slight variation at the beginning due to different file names.   

```
bl-filter-illumina -a -i ../../00-DATA/fastq/SRR[number]_pass_1.fastq.gz -i ../../00-DATA/fastq/[number]_pass_2.fastq.gz -o SRR[number].1.fq -o SRR[number].2.fq -u SRR[number].unp.fq > blf.out 2> blf.err &
```

Copy the trimmed adapter file and concatonate the unpaired sequences to a new file. Link the second trimmed adapter file to a new file so that file naming was consistent. 

```
cp SRR[number].1.fq all_1.fq
```

```
cat SRR[number].unp.fq >> all_1.fq
```

```
ln -s SRR2830762.2.fq.gz all_2.fq.gz
```

2.2.2 Use the script ```fix_names.pl``` in order to fix the deflines of the SRA files so that they are formatted properly for transcriptome assembly with Trinity. This script is available in the scripts directory in this repository. 

```
perl /bwdata1/jfryan/38-SIMION_DATA/02-FILTER_ILLUMINA/fix_names.pl SRR2484238.1.fq.gz SRR2484238.2.fq.gz SRR2484238.unp.fq.gz > fix_sra_names.out 2> fix_sra_names.err
```

#### 2.3 De novo transcriptome assembly with Trinity v2.4.0.  

```
/usr/local/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --max_memory 750G --CPU 12 --left ../01-BL-FILTER/SRR[number].1.fq.renamed --right ../01-BL-FILTER/SRR[number].2.fq.renamed --full_cleanup --normalize_reads --normalize_max_read_cov 30 > trin.out 2> trin.err &
```
##### Note 12/21/17: Trinity runs were initially started with higher memory and central processors, however, we lowered those values due to server constraints.

#### 2.4 Use RSEM (Li and Dewey, 2011) to measure the gene and isoform abundance.
```
/usr/local/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts ../02-TRINITY/trinity_out_dir.Trinity.fasta --seqType fq --left ../01-BL-FILTER/SRR5755244.1.fq.gz --right ../01-BL-FILTER/SRR5755244.2.fq.gz --output_dir aea --est_method RSEM --aln_method bowtie2 --thread_count 100 --prep_reference > aea.out 2> aea.err &
```

```
rsemgetbestseqs.py ./aea/RSEM.isoforms.results ../02-TRINITY/trinity_out_dir.Trinity.fasta > rgbs.out 2> rgbs.err &
```

#### 2.5 Translate holothurian nucleotide transcriptome sequences into amino acid sequences with TransDecoder v3.0.0. We will set the –m flag to 50 and use the results from blast and hmmscan searches to inform the final TransDecoder prediction step.  

```
TransDecoder.LongOrfs -t [transcriptome_file] -m 50 > td.out 2> td.err
```

```
blastp -query longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > blastp.out 2> blastp.err 
```

```
hmmscan --cpu 1 --domtblout outfile.domtblout Pfam-A.hmm longest_orfs.pep > hs.out 2> hs.err
```

```
TransDecoder.Predict -t [transcriptome_file] --retain_pfam_hits out.domtblout --retain_blastp_hits out.blastp.out > tdp.out 2> tdp.err
```

#### 2.6 The program [Alien Index](https://github.com/josephryan/alien_index) and a database of representative metazoan and non-metazoan sequences (http://ryanlab.whitney.ufl.edu/downloads/alien_index/) will allow us to remove any contaminating, non-metazoan sequences. 

```
blastp -query [infile.pep.fa] -db ai.fa -outfmt 6 -max_target_seqs 1000 -seg yes -evalue 0.001 -out [file.out] > file.std 2> file.err
```

```
../alien_index --blast=[file_ai.out] --alien_pattern=ALIEN [out.alien_index] > ai.out 2> ai.err 
```

```
perl /usr/local/bin/remove_aliens.pl [out.alien_index] [original_transcriptome.fa] > [filtered_transcriptome.fa] > ra.out 2> ra.err
```

#### 2.7 Identify orthogroups across holothurian transcriptomes with OrthoFinder v1.1.8. 
```
orthofinder -f [dir_w_protein_fast_files] -op > of.out 2> of.err
```
##### Note 12/21/17: The first step generated 676 blastp files. We used the script `blastp_parser.pl` [which can be found in the script repository] in order to separate the blastp files into smaller .sh files so they could be executed sequentially on the servers. 

```
blastp -outfmt 6 -evalue 0.001 -query [renamed_fasta_file_w_all_seqs] -db BlastDB -out outfile.txt > blastp.out 2> blastp.err
```

```
orthofinder -b [dir_w_blast_results] > ofb.out 2> ofb.err
```
##### Note: 12/21/17: At this point we recognized a mis-labeled transcriptome from the NCBI database indicating that there was a duplicate transcriptome. I re-ran: 
```
orthofinder -b [dir_w_blast_results] > ofb.out 2> ofb.err
```
##### and removed the duplicate transcriptome by placing a `#` in front of the duplicate species in the `SpeciesIDs.txt` file generated from the previous analysis. 

```
python trees_from_MSA.py [dir_w_orthofinder_results] > tfm.out 2> tfm.err
```

#### 2.8  

#### 2.4 Generate single copy orthogroups. First, the script ```filter_ogs_write_scripts.pl``` (available in the scripts directory of this repository) retains orthogroup fasta files that contain a user-specified minimum number of taxa (need to determine the scope of this). Lastly, ```filter_ogs_write_scripts.pl``` automates the following processes: 

2.4.1 sequences within each orthogroup are aligned using Mafft v7.309 

```mafft-linsi --localpair --maxiterate 1000 --thread 20 [infile] >mafft.out 2> mafft.err```

2.4.2 alignments are refined using Gblockswrapper v0.03 (https://goo.gl/fDjan6)

```Gblockswrapper [infile.mafft] > outfile.mafft-gb > gbw.out 2> gbw.err```

2.4.3 Gblockswrapper sometimes leaves blank sequences that cause downstream issues; the ```remove_empty_seqs``` script, available in the scripts directory of this repository, removes empty sequences and spaces from sequence lines. 

```remove_empty_seqs [outfile.mafft-gb] > res.out 2> res.err```

2.4.4 Maximum-likelihood orthogroup gene trees are estimated in IQTree v1.5.5 

```iqtree-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m LG -pre [output prefix] > iq.out 2> iq.err```

2.4.5 orthogroups with multiple *M. ovum* sequences are pruned in PhyloTreePruner v1.0 

```java PhyloTreePruner [infile.tree] 28 [infile.align] 0.5 u > ptp.out 2> ptp.err```


#### 2.5 Concatenate (insert #) single-copy loci filtered from step 5 to create a matrix and partition file for use in downstream phylogenomic analyses using ```fasta2phylomatrix``` (available in the scripts directory of this repository). Definition lines in each fasta file were edited (```perl -pi.orig -e 's/\|.*$//;' *.fa```) prior to running ```fasta2phylomatrix```.  

#### 2.6 Estimate species phylogeny using concatenated and coalescent gene tree/species tree methods.  

2.6.1 Concatenated matrix, Maximum Likelihood: estimate a bootstrapped (1000 ultrafast replicates) species phylogeny in IQtree v1.5.5 using the concatenated dataset. We will use the flag -m TEST to find best partition scheme and estimate the tree. The partition file will be created with the script ```fasta2phylomatrix```, which is available in this respository.

```
iqtree-omp –nt [#threads] –s [infile] –pre [outfile_prefix] –spp [partition file] –m TEST –bb 1000 –bspec GENESITE > iqo.out 2> iqo.err
```

2.6.2 Concatenated matrix, Bayesian inference: estimate species phylogeny in PhyloBayes-MPI v1.7 using the concatenated dataset. If PhyloBayes is not close to convergence after 1 month runtime, we will use the jackknife approach described in Simion et al. 2017.  

```
mpirun -n [# threads] pb_mpi -d [infile.phy] -cat -gtr chain1 > pb1.out 2> pb1.err
mpirun -n [# threads] pb_mpi -d [infile.phy] -cat -gtr chain2 > pb2.out 2> pb2.err
bpcomp -x [burnin] [sample_every_x_number_of_trees] <chain1> <chain2> > bpcomp.out 2> bpcomp.err
```

2.6.3 Coalescent-based phylogeny: estimate the species phylogeny using ASTRAL-II v4.11.1 and ASTRID v1.4. 

> i) Generate individual maximum-likelihood gene trees in IQtree. 

```
iqtree-omp –nt AUTO –s [infile] –pre [prefix_for_outfiles] –m MFP+MERGE –bb 1000 > iq.out 2> iq.err
```

> ii) ASTRAL-II constrains the search space to those species trees that derive their bipartitions from the input gene trees

```
java -Xmx1000000M -jar astral.jar -i [gene_trees_file] -o [output_file] > astral.out 2> astral.err
```

> iii) ASTRID uses a distance matrix generated from the input gene trees to estimate the species tree and is robust to missing data

```
ASTRID –i [infile] –o [outfile] –m bionj > astrid.out 2> astrid.err
```

> iv) Compute branch support using local posterior probabilities.  

#### 2.7 If there are conflicting species-tree topologies from 2.7, perform SOWH tests (implemented in sowhat v.0.36) to compare topologies. Any topologies that can be rejected with a P-Value <= 0.05 will be excluded from downstream analyes (but still reported in results). Constraint trees will be added to the phylotocol before running the tests.

2.7.1 example sowhat command line

```sowhat --constraint=[topology_to_be_tested] --aln=[alignment] --name=[name] --dir=[output_dir] --rax=[raxmlHPC-PTHREADS-SSE3 -T [num_threads]] ```


## 3 WORK COMPLETED SO FAR WITH DATES  

October 5 2017- we have completed steps [insert] prior to release of phylotocol version 1.0
