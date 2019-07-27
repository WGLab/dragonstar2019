## Before the tutorial

1. To use installed softwares, you need to run 
```
. /shared/miniconda3/etc/profile.d/conda.sh
conda activate base
```
2. go to your home folder, and `mkdir -p ~/project/alignment/` folder to create an alignment folder

## The tutorial for read alignment and variants calling

### 1. Short read alignment and variants calling
#### 1.1 Preparation of the folder and data
1. `cd ~/project/alignment/`, and then `mkdir short-reads` and `cd short-reads`
2. Link data by `ln -s /shared/data/NA12878_test_data/short-reads data`
3. Link reference genome by `ln -s /shared/data/ref_hg37_chr1/ref ref`
4. Link high quality variant file by `ln -s /shared/data/ref_hg37_chr1/vcf vcf`

#### 1.2 Short reads alignment
In this tutorial, the reads are from a 2MB region in chr1. We will test two ways to do the alignment.

Fist, we check what input data we will use:

```
(base) [biouser14@main-lx short-reads]$ wc -l data/sr.chr1.2mb_1.fq data/sr.chr1.2mb_2.fq
  2432672 data/sr.chr1.2mb_1.fq
  2432672 data/sr.chr1.2mb_2.fq
  4865344 total
```

So there are about 1.2 million paired-end reads (remember that each read uses 4 lines in the FASTQ file). 

You can also check the read length yourself, by examining the first a few lines of the FASTQ file. You will see that it is 148bp in each read.

##### 1.2.1 Alignment with minimap2. 
```
minimap2 -ax sr ref/hg37d5.chr1.fa data/sr.chr1.2mb_1.fq data/sr.chr1.2mb_2.fq | samtools sort | samtools view -bS - > chr1.2mb.mp2.bam
samtools index chr1.2mb.mp2.bam
```
Reads will be aligned with chr1, and then sorted and saved into a bam file and then build index for it.

The first command should finish in a few minutes, and should use less than 4GB memory. The second command should finish instantly. If your command takes very long time to run, please report to TA to solve the issue.

You can read the manual for minimap2 [here](https://lh3.github.io/minimap2/minimap2.html). As you will see, we used `-ax sr` for "short genomic paired-end reads" because this is a short-read sequencing data set with 148bp read length. When two read files are specified, minimap2 reads from each file in turn and merge them into an interleaved stream internally. Two reads are considered to be paired if they are adjacent in the input stream and have the same name.

The second command is used to build an index, and you can see that a new file `chr1.2mb.mp2.bam.bai` is generated. The index file is important to facilitate downstream analysis on the BAM file.

Note that here we used the `hg37d5.chr1.fa` file as reference, which contains sequence for chromosome 1 only. This allows the software to complete the alignment in a short period of time. In real-world settings, you should align reads to the entire genome.


##### 1.2.2 Alignment with bwa-mem. 
```
bwa mem ref/hg37d5.chr1.fa data/sr.chr1.2mb_1.fq data/sr.chr1.2mb_2.fq | samtools sort | samtools view -bS - > chr1.2mb.bwa.bam
samtools index chr1.2mb.bwa.bam
```
An alignment bam file will be generated and indexed for further analysis.

Similar to the minimap2 command used in the procedure above, here we only align to `hg37d5.chr1.fa` (chromosome 1), to speed up the exercise. However, you will see that this command takes much longer time to finish, compared to minimap2. Generally speaking, it should finish within 5-10 minutes. You can check the screen for progress.

You can read the manual for bwa-mem [here](http://bio-bwa.sourceforge.net/bwa.shtml). As you will see, we supplied two FASTQ files in the command, so this command assumes the i-th read in the first file `sr.chr1.2mb_1.fq` and the i-th read in the second file `sr.chr1.2mb_2.fq` constitute a read pair. In the paired-end mode, the mem command will infer the read orientation and the insert size distribution from a batch of reads.

Since the cloud server has two available CPU cores per node, we can speed up the alignment by using multiple threads. You can try the command below to see how it speeds up the alignment:

```
bwa mem -t 2 ref/hg37d5.chr1.fa data/sr.chr1.2mb_1.fq data/sr.chr1.2mb_2.fq | samtools sort | samtools view -bS - > chr1.2mb.bwa.bam
```

When running the program, you can open a second terminal window, and run the `top` command to see the CPU usage (the display refreshes every a few seconds, and you can press "q" to quite the `top` command). Since two cores are requested, you will see that the CPU usage (%CPU) is 200. In practice, most whole genome sequencing studies nowadays use multiple cores (and likely multiple computing nodes) to speed up the alignment procedure.

#### 1.3 View bam files

We can examine the output bam file in a few different ways.

##### 1.3.1 View bam records
```
samtools view chr1.2mb.bwa.bam | less
``` 
to see how to represent alignment for each read. Detailed description of the SAM format is described [here](https://samtools.github.io/hts-specs/SAMv1.pdf). You can type `q` to exit the `less` session.

Alternatively, you can use `samtools view chr1.2mb.bwa.bam | head` to print the first 10 lines to the terminal, and then examine the output yourself.

###### 1.3.2 tview bam files
```
samtools tview -p 1:156000000 chr1.2mb.bwa.bam ref/hg37d5.chr1.fa
``` 
to see how alignments parallel with the reference genome (the first line in the screen).

You can press left and right key in the keyboard to move to the left or right of the screen. 

To jump to a specific genomic position, you can type `g`, and then enter the coordiante. Try to jump to the "1:156000211" position.

To quite the tview session, just press `q`.


#### 1.4 Variants calling
A simple way for variant calling is to use bcftools. We can generate variant calling for both bam files generated above.

1. On bam file aligned with `minimap2`. 
```
bcftools mpileup -r 1:155000000-157000000 -f ref/hg37d5.chr1.fa chr1.2mb.mp2.bam | bcftools call -mv -Ob -o mp2.bcftools.call.bcf
```
will call variants and save to `mp2.bcftools.call.bcf` in a bcf format. 

You can view the file using 
```
bcftools view mp2.bcftools.call.bcf | less
``` 
or 
```
bcftools view -i '%QUAL>=20' mp2.bcftools.call.bcf | less
``` 
to only view those higher quality (>20) variants.

BCF is a binary format, and sometimes you may want to have a text-format VCF file for additional analysis and examination. If you want to convert the file formats, you can use

```
bcftools view mp2.bcftools.call.bcf > mp2.bcftools.call.vcf
``` 
to convert bcf to vcf.

2. On A bam file aligned with `bwa mem`
Simiarly, 
```
bcftools mpileup -r 1:155000000-157000000 -f ref/hg37d5.chr1.fa chr1.2mb.bwa.bam | bcftools call -mv -Ob -o bwa.bcftools.call.bcf
``` 
to call variants and save to `bwa.bcftools.call.bcf`. 

The first mpileup part generates genotype likelihoods at each genomic position with coverage, and use a pipe to feed the output to the second part; the second call part makes the actual calls. You can add `--threads 2` to the second command to use two threads, but since the more time consuming step is the first step, it will not shorten the waiting time much.

Generally speaking, the command should finish within 5-10 minutes. Note that we restricted pileup to a specific genomic region, and used a reference file for chromosome 1 only, so that we can shorten the waiting time for this command.

You can also use `bcf view` to what variants have been called.

Also, 
```
bcftools view bwa.bcftools.call.bcf > bwa.bcftools.call.vcf
``` 
can be used to convert bcf to vcf for further analysis.

#### 1.5 Evaluation of variants calling

We next generate statistics of the called variants. Note that we used default parameters in bcftools for variant calling, so the variant calls are not optimal. More specialized variant calling tools (such as GATK and Strelka2) that use population information or allow haplotype-based calling are available and can be used to improve accuracy.

##### 1.5.1. On bam file aligned with `minimap2`. 
```
bcftools view -i '%QUAL>=200' mp2.bcftools.call.bcf > mp2.bcftools.call.vcf
bgzip -f mp2.bcftools.call.vcf -c > mp2.bcftools.call.vcf.gz
bcftools index -f mp2.bcftools.call.vcf.gz
bcftools stats mp2.bcftools.call.vcf.gz | grep "^SN"
```
will generate the statistics below
```
SN      0       number of samples:      1
SN      0       number of records:      2060
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1765
SN      0       number of MNPs: 0
SN      0       number of indels:       295
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```

In the commands above, we first filter a subset of more confident variant calls and write them to a VCF file `mp2.bcftools.call.vcf`, then compress the file and perform some simple statistics on the file using `bcftools`.

We can see that there are 1765 snps and 295 indels among 2 millions bp region. Generally speaking, we expect to see 1 variant per 1000 base pairs in human genome, so this is consistent with our expectation.

```
bcftools stats mp2.bcftools.call.vcf.gz | grep "TSTV"
```
will give the statistics below:
```
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       1260    505     2.50    1260    505     2.50
```
where TS is 2.5 times compared with TV. According to [wiki](https://en.wikipedia.org/wiki/Transition_%28genetics%29), "Transition, in genetics and molecular biology, refers to a point mutation that changes a purine nucleotide to another purine (A ↔ G), or a pyrimidine nucleotide to another pyrimidine (C ↔ T)". Generall speaking, for a whole-genome sequencing data set, the Transition/Transversion ratio should be around 2, and for a whole-exome sequencnig data set, the ratio should be around 3. Here we are looking at a very small genomic region with a small number of variants, so the statistics may deviate from whole genome expectation.

###### 1.5.1.1 Evaluation of the variant calling

We next compare the called variants against high-quality variants in a gold-standard set.
```
bcftools isec mp2.bcftools.call.vcf.gz vcf/vcf.chr1.2mb.vcf.gz -p minimap2_perf
```

The gold-standard set consists of high quality variant calls retrieved from the [genome in a bottle consortium](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/) on this particular genomic region.

The first command will generate the intersection of called variants and the gold-standard variants, and you have a new folder called `minimap2_perf` with the files below
```
minimap2_perf/0000.vcf  for records private to  mp2.bcftools.call.vcf.gz
minimap2_perf/0001.vcf  for records private to  vcf/vcf.chr1.2mb.vcf.gz
minimap2_perf/0002.vcf  for records from mp2.bcftools.call.vcf.gz shared by both        mp2.bcftools.call.vcf.gz vcf/vcf.chr1.2mb.vcf.gz
minimap2_perf/0003.vcf  for records from vcf/vcf.chr1.2mb.vcf.gz shared by both mp2.bcftools.call.vcf.gz vcf/vcf.chr1.2mb.vcf.gz
```

We thus investigate `minimap2_perf/0002.vcf` to see the intersected variants.
```
bcftools stats minimap2_perf/0002.vcf | grep "^SN"
```
And then,
```
SN      0       number of samples:      1
SN      0       number of records:      1588
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1569
SN      0       number of MNPs: 0
SN      0       number of indels:       19
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
It seems that there are 1588 called variants which are consistent with gold standard, and there are total 2060 called variants, and thus the precision is 1588/2060=0.771. If only SNPs are considered, the precision is 1569/1765=0.89.

Using 
```
bcftools stats vcf/vcf.chr1.2mb.vcf.gz | grep "^SN"
```
One will have 
```
SN      0       number of samples:      1
SN      0       number of records:      1922
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1587
SN      0       number of MNPs: 0
SN      0       number of indels:       336
SN      0       number of others:       0
SN      0       number of multiallelic sites:   27
SN      0       number of multiallelic SNP sites:       1
```
We can know that there are 1922 gold-standard variants and among them, 1587 are SNPs, and thus, the recall is 1588/1922=0.826, and the recall of SNPs is 1569/1587=0.989. Note that there are also a few multiallelic sites (sites with two or more alternative alleles) that we did not include in the calculation.


##### 1.5.2. On A bam file aligned with `bwa bam`
```
bcftools view -i '%QUAL>=200' bwa.bcftools.call.bcf > bwa.bcftools.call.vcf
bgzip -f bwa.bcftools.call.vcf -c > bwa.bcftools.call.vcf.gz
bcftools index -f bwa.bcftools.call.vcf.gz
bcftools stats bwa.bcftools.call.vcf.gz | grep "^SN"
```
will generate the statistics below
```
SN      0       number of samples:      1
SN      0       number of records:      2074
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1772
SN      0       number of MNPs: 0
SN      0       number of indels:       302
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
where there are 1772 snps and 302 indels among 2 millions bp region. 

```
bcftools stats bwa.bcftools.call.vcf.gz | grep "TSTV"
```
will give the statistics below:
```
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       1266    506     2.50    1266    506     2.50
```
where TS is 2.5 times compared with TV.

###### 1.5.2.1 Evaluation of the variant calling
This called variants is also compared against the gold-standard variants, and precision and recall will be calculated.
```
bcftools isec bwa.bcftools.call.vcf.gz vcf/vcf.chr1.2mb.vcf.gz -p bwa_perf
bcftools stats bwa_perf/0002.vcf | grep "^SN"
```
Similarly, we get
```
SN      0       number of samples:      1
SN      0       number of records:      1590
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1571
SN      0       number of MNPs: 0
SN      0       number of indels:       19
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
The precision is 1590/2074=0.767, and the recall is 1590/1922=0.827. If only SNPs are considered, the precision is 1571/1772=0.887, and the recall is 1571/1587=0.99.


#### 1.6 Comparison of two called variants.
We next to see how the variants overlap between the two called variants by the two tools.
```
bcftools isec mp2.bcftools.call.vcf.gz bwa.bcftools.call.vcf.gz -p twoshort_comparison
bcftools stats twoshort_comparison/0002.vcf | grep "^SN"
```

and you will have 
```
SN      0       number of samples:      1
SN      0       number of records:      2041
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1764
SN      0       number of MNPs: 0
SN      0       number of indels:       277
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
Where you can see that the vast majority of the two sets of called variants are same. So in this case, variant callers have much higher influence on the final results than the aligners.

We should emphasize that the above comparative analysis is a quick-and-dirty way of checking variant overlaps using `bcftools`. We did not consider many other issues, such as indel re-alignment, the splitting of multi-allelic variants, the heterozygosity of the variants, the regions with low mappability, and so on. For more formal comparison of varint call sets, there are more specialized tools, such as the [vcfeval](https://github.com/RealTimeGenomics/rtg-tools) and [hap.py/som.py](https://github.com/Illumina/hap.py) released by Illumina. For a description on how GA4GH evaluates variant calling performance using GIAB data, you can check [here](https://github.com/ga4gh/benchmarking-tools/blob/master/resources/high-confidence-sets/giab.md).

### 2. Long read alignment and variants calling

#### 2.1 Preparation of the folder and data
1. `cd ~/project/alignment/`, and then `mkdir long-reads` and `cd long-reads`
2. Link data by `ln -s /shared/data/NA12878_test_data/long-reads data`
3. Link reference genome by `ln -s /shared/data/ref_hg37_chr1/ref ref`
4. Link high quality variant file by `ln -s /shared/data/ref_hg37_chr1/vcf vcf`

#### 2.2 Long reads alignment
```
minimap2 --MD -ax map-ont ref/hg37d5.chr1.fa data/np.chr1.2mb.fq  | samtools sort | samtools view -bS - > chr1.2mb.bam
samtools index chr1.2mb.bam
```
The commands above will align long reads in `data/np.chr1.2mb.fq`, and then sort/save alignment into a bam file `chr1.2mb.bam`. 
An index is also created so that you can use `samtools view` or `samtools tview` to view the alignment in the bam file.


#### 2.3 View bam files
One can try one of the commands below to view the bam files.
```
samtools view chr1.2mb.bam | less
samtools tview -p 1:156000000 chr1.2mb.bam ref/hg37d5.chr1.fa
```

#### 2.4 Variants calling
A simple tool `longshot` can be used to call variants from long-reads aligned bam file. To use this tool, you need to `conda activate longshot` to activate the virtual environment.

Then,
```
longshot --bam chr1.2mb.bam --ref ref/hg37d5.chr1.fa --out mp2.longshot.o.vcf
``` 
will generate called variants and save in `mp2.longshot.o.vcf` in a vcf format. 

VCF format is plain text, and you can use `less mp2.longshot.o.vcf` to see what is inside this file.

#### 2.5 Statistics of called variants
We also want to see how many snps and indels are in the vcf files.
```
bcftools view -i '%QUAL>=30' mp2.longshot.o.vcf > mp2.longshot.vcf
bgzip -f mp2.longshot.vcf -c > mp2.longshot.vcf.gz
bcftools index -f mp2.longshot.vcf.gz
bcftools stats mp2.longshot.vcf.gz | grep "TSTV"
```
will tell you 
```
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       1616    529     3.05    1616    529     3.05
```
Where TS is more than 3 times than TV.

```
bcftools stats mp2.longshot.vcf.gz | grep "^SN"
```
will tell that 
```
SN      0       number of samples:      1
SN      0       number of records:      2145
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 2145
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
Where there are 2145 snps and no indels found. 

#### 2.5.1. Evaluation of called variants.
We also validate the called variants with the gold-standard variants. 
```
bcftools isec mp2.longshot.vcf.gz vcf/vcf.chr1.2mb.vcf.gz -p longshot_perf
bcftools stats longshot_perf/0002.vcf | grep "^SN"
```
And then, we can know that 
```
SN      0       number of samples:      1
SN      0       number of records:      1189
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1189
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
Since no indels are detected by `longshot`, we only calculate recall and precision for snps: the recall is 1189/1587=0.75, but the precision is 1189/2145=0.55, which is much worse than the performance generated by short reads. 


#### 2.6 Comparison of variants called from short-reads and variants called from long reads
Since there is a big difference in called variants from long reads compared with variants called from short reads, we will also investigate the insection between the variants called from short reads and long reads using the command below:
```
bcftools isec ../short-reads/mp2.bcftools.call.vcf.gz mp2.longshot.vcf.gz -p shortlong_comparison
```

Then, we check the overlapped variants using
```
bcftools stats shortlong_comparison/0002.vcf | grep "^SN"
```
where you can see 
```
SN      0       number of samples:      1
SN      0       number of records:      1294
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1294
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
There are 1294 intersected variants called from short reads and long reads. 

In general, due to the high error rate in base calling, more sophisticated variant calling tools would be needed to address challenges in long-read data for SNP/indel calling, such as [DeepVariant](https://github.com/google/deepvariant) developed by Google and [Clairvoyante](https://github.com/aquaskyline/Clairvoyante) developed by Ruibang Luo.


## After the tutorial

To do other tutorial, you might need to run `conda deactivate` to go back to the base environment for other projects. ***If you still have issues to run other projects, please re-login.***

