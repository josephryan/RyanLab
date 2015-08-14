## Allpaths-LG

<code>ErrorCorrectReads.pl PAIRED_READS_A_IN=Ef.1.fq PAIRED_READS_B_IN=Ef.2.fq PAIRED_SEP=100 THREADS=10 PHRED_ENCODING=64 READS_OUT=Efet</code>

## Cutadapt

<code>cutadapt -f fastq -O 6 -q 20 -a AGATCGGAAGAGC -o Efet.A.fq ../01-Clean/Efet.paired.A.fastq</code>

## SOAPdenovo

<code>SOAPdenovo-63mer all -K31 -p 5 -o Efet.Soap31 -s config</code>

<pre>#config file
max_rd_len=100
[LIB]
avg_ins=300
reverse_seq=0
asm_flags=3
rd_len_cutoff=100
rank=1
pair_num_cutoff=3
map_len=32
q1=../02-Cutadapt/Efet.A.fq
q2=../02-Cutadapt/Efet.B.fq
q=../02-Cutadapt/Efet.unp.fq</pre>

## ABySS

<code>abyss-pe j=5 k=31 name=Efet31 in='Efet.A.fq Efet.B.fq' se=Efet.unp.fq</code>

## Platanus

<code>platanus assemble -o Efetplat32 -k 32 -m 500 -f Efet.A.fq Efet.B.fq Efet.unp.fq</code>

<code>platanus scaffold -o Efetplat32.scaf -c Efetplat32_contig.fa -b Efetplat32_contigBubble.fa -IP1 Efet.A.fq Efet.B.fq</code>

<code>platanus gap_close -c Efetplat32.scaf_scaffold.fa -IP1 Efet.A.fq Efet.B.fq</code>

## CEGMA

<code>cegma –genome assembly.fa kogs.fa</code>

## BLAT

<code>blat Efet.Soap31.scafSeq Efet_EST.fa Soap31.blat</code>

## Isoblat

<code>isoblat Soap31.blat Efet_EST.fa</code>

## HMMsearch

<code>hmmsearch -A EfetHD.aln -o EfetHD.hmmsearch hd60.hmm Efet.aa</code>

## MAFFT

<code>mafft --globalpair --maxiterate 1000 --clustalout Efet_sequences.phy > Efet_alignment2.phy</code>

## RAxML

<code>raxmlHPC-PTHREADS -T 30 -s Spiralia.phy -n Spiralia -m PROTGAMMALG -p 1234</code>

## MrBayes

<code>mb file.nex</code>

# file.nex included an alignment and this code block
<pre>BEGIN MRBAYES;
prset aamodelpr=mixed; lset rates = gamma; mcmcp mcmcdiagn = no nruns = 2 ngen = 5000000 printfreq = 5000 samplefreq = 500 nchains = 5 savebrlens = yes; mcmc; sumt filename =  nRuns = 1 Relburnin = YES BurninFrac = .25 Contype = Allcompat;).
END;</pre>



## LINKS TO SOFTWARE PROGRAMS

* Allpaths - https://www.broadinstitute.org/software/allpaths-lg/blog/?page_id=12
* Cutadapt - https://code.google.com/p/cutadapt/
* SOAPdenovo - http://soap.genomics.org.cn/soapdenovo.html
* ABySS - http://www.bcgsc.ca/platform/bioinfo/software/abyss
* Platanus - http://platanus.bio.titech.ac.jp
* CEGMA - http://korflab.ucdavis.edu/datasets/cegma/
* BLAT - http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads
* Isoblat - https://github.com/josephryan/isoblat
* MAFFT - http://mafft.cbrc.jp/alignment/software/source.html
* RAxML - http://sco.h-its.org/exelixis/web/software/raxml/i
* MrBayes - http://mrbayes.sourceforge.net/download.php
