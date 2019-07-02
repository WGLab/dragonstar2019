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
In this tutorial, the reads are from a 2MB region in chr1. There are two ways to do the alignment.

##### 1.2.1 Alignment with minimap2. 
```
minimap2 -ax sr ref/hg37d5.chr1.fa data/sr.chr1.2mb_1.fq data/sr.chr1.2mb_2.fq | samtools sort | samtools view -bS - > chr1.2mb.mp2.bam
samtools index chr1.2mb.mp2.bam
```
Reads will be aligned with chr1, and then sorted and saved into a bam file and then build index for it.

##### 1.2.2 Alignment with bwa-mem. 
```
bwa mem ref/hg37d5.chr1.fa data/sr.chr1.2mb_1.fq data/sr.chr1.2mb_2.fq | samtools sort | samtools view -bS - > chr1.2mb.bwa.bam
samtools index chr1.2mb.bwa.bam
```
Ann alignment bam file will be generated and indexed for further analysis.

#### 1.3 View bam files
One can see what is bam in two ways.

##### 1.3.1 View bam records
```
samtools view chr1.2mb.bwa.bam | less
``` 
to see how to represent alignment for each reads. 
`type q` to exit.

###### 1.3.2 tview bam files
```
samtools tview -p chr1:156000000 chr1.2mb.bwa.bam ref/hg37d5.chr1.fa
``` 
to see how alignments parallel with the reference genome.

#### 1.4 Variants calling
A simple way for variant calling is to use bcftools. We can generate variant calling for both bam files generatd above.

1. On bam file aligned with `minimap2`. 
```
bcftools mpileup -f ref/hg37d5.chr1.fa chr1.2mb.mp2.bam | bcftools call -mv -Ob -o mp2.bcftools.call.bcf
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

If you want to need vcf files for analysis, please use
```
bcftools view mp2.bcftools.call.bcf > mp2.bcftools.call.vcf
``` 
to convert bcf to vcf;

2. On A bam file aligned with `bwa bam`
Simiarly, 
```
bcftools mpileup -f ref/hg37d5.chr1.fa chr1.2mb.bwa.bam | bcftools call -mv -Ob -o bwa.bcftools.call.bcf
``` 
to call variants and save to `bwa.bcftools.call.bcf`. 
You can also use `bcf view` to what variants have been called.

Also, 
```
bcftools view bwa.bcftools.call.bcf > bwa.bcftools.call.vcf
``` 
can be used to convert bcf to vcf for further analysis.

#### 1.5 Statistics and evaluation of variants calling
We next generate statistics of the called variants.

##### 1.5.1. On bam file aligned with `minimap2`. 
```
bgzip mp2.bcftools.call.vcf
bcftools index -f mp2.bcftools.call.vcf.gz
bcftools stats mp2.bcftools.call.vcf.gz | grep "^SN"
```
will generate the statistics below
```
SN      0       number of samples:      1
SN      0       number of records:      1248
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1103
SN      0       number of MNPs: 0
SN      0       number of indels:       145
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
where there are 1103 snps and 145 indels among 1 millions bp region. The indels are less than 14 percentage of snps. 

```
bcftools stats mp2.bcftools.call.vcf.gz | grep "TSTV"
```
will give the statistics below:
```
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       762     341     2.23    762     341     2.23
```
where TS is 2.23 times compared with TV. According to [wiki](https://en.wikipedia.org/wiki/Transition_%28genetics%29), "Transition, in genetics and molecular biology, refers to a point mutation that changes a purine nucleotide to another purine (A ↔ G), or a pyrimidine nucleotide to another pyrimidine (C ↔ T)"

###### 1.5.1.1 Evaluation of the variant calling
We next compare the called variants against high-quality variants in a gold-standard set.
```
bcftools isec mp2.bcftools.call.vcf.gz vcf/vcf.chr1.2mb.vcf.gz -p minimap2_perf
```
The first command will generate the intersection of called variants and the gold-standard variants, and you have a new folder called `minimap2_perf` with the files below
```
minimap2_perf/0000.vcf    for records private to  mp2.bcftools.call.vcf.gz
minimap2_perf/0001.vcf    for records private to  bwa.bcftools.call.vcf.gz
minimap2_perf/0002.vcf    for records from mp2.bcftools.call.vcf.gz shared by both        mp2.bcftools.call.vcf.gz bwa.bcftools.call.vcf.gz
minimap2_perf/0003.vcf    for records from bwa.bcftools.call.vcf.gz shared by both        mp2.bcftools.call.vcf.gz bwa.bcftools.call.vcf.gz
```

We thus investigate `minimap2_perf/0002.vcf` to see the intersected variants.
```
bcftools stats minimap2_perf/0002.vcf | grep "^SN"
```
And then,
```
SN      0       number of samples:      1
SN      0       number of records:      234
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 233
SN      0       number of MNPs: 0
SN      0       number of indels:       1
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
It seems that there are 234 called variants which are correct, and there are total 1248 called variants, and thus the precision is 234/1248=0.1875. Using 
```
bcftools stats vcf/vcf.chr1.2mb.vcf.gz | grep "^SN"
```
We can know that there are 829 gold-standard variants, and thus, the recall is 234/829=0.28.


##### 1.5.2. On A bam file aligned with `bwa bam`
```
bgzip bwa.bcftools.call.vcf
bcftools index -f bwa.bcftools.call.vcf.gz
bcftools stats bwa.bcftools.call.vcf.gz | grep "^SN"
```
will generate the statistics below
```
SN      0       number of samples:      1
SN      0       number of records:      1415
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1249
SN      0       number of MNPs: 0
SN      0       number of indels:       166
SN      0       number of others:       0
SN      0       number of multiallelic sites:   2
SN      0       number of multiallelic SNP sites:       0
```
where there are 1249 snps and 166 indels among 1 millions bp region. The indels are less than 15 percentage of snps. 

```
bcftools stats bwa.bcftools.call.vcf.gz | grep "TSTV"
```
will give the statistics below:
```
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       847     402     2.11    847     402     2.11
```
where TS is 2.11 times compared with TV.

###### 1.5.2.1 Evaluation of the variant calling
This called variants is also compared against the gold-standard variants, and precision and recall will be calculated.
``
bcftools isec bwa.bcftools.call.vcf.gz vcf/vcf.chr1.2mb.vcf.gz -p bwa_perf
bcftools stats bwa_perf/0002.vcf | grep "^SN"
``
Similarly, we get
```
SN      0       number of samples:      1
SN      0       number of records:      362
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 360
SN      0       number of MNPs: 0
SN      0       number of indels:       2
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
The precision is 362/1415=0.256, and the recall is 362/829=0.437. Compared with the variants called with `minimap2`, it seems that on short reads, `bwa` gives better performance.


#### 1.6 Comparison of two called variants.
We next to see how the variants overlap between the two called variants by the two tools.
```
bcftools isec mp2.bcftools.call.vcf.gz bwa.bcftools.call.vcf.gz -p twoshort_comparison
bcftools stats twoshort_comparison/0002.vcf | grep "^SN"
```

and you will have 
```
SN      0       number of samples:      1
SN      0       number of records:      901
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 791
SN      0       number of MNPs: 0
SN      0       number of indels:       110
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
Where you can see that about 30% variants are different between the two sets of called variants. But one can try to use quality filer `-i '%QUAL>=20'` to investigate high-quality varriants only to see whether higher overlapping percentage can be obtained.


### 2. Long read alignment and variants calling
#### 2.1 Preparation of the folder and data
1. `cd ~/project/alignment/`, and then `mkdir long-reads` and `cd long-reads`
2. Link data by `ln -s /shared/data/NA12878_test_data/long-reads data`
3. Link reference genome by `ln -s /shared/data/ref_hg37_chr1/ref ref`
4. Link high quality variant file by `ln -s /shared/data/ref_hg37_chr1/vcf vcf`

#### 2.2 Long reads alignment
```
minimap2 -ax map-ont ref/hg37d5.chr1.fa data/np.chr1.2mb.fq  | samtools sort | samtools view -bS - > chr1.2mb.bam
samtools index chr1.2mb.bam
```
The commands above will align long reads in `data/np.chr1.2mb.fq`, and then sort/save alignment into a bam file `chr1.2mb.bam`. 
An index is also created so that you can use `samtools view` or `samtools tview` to view the alignment in the bam file.

#### 2.3 View bam files
One can try one of the commands below to view the bam files.
```
samtools view chr1.2mb.bam | less
samtools tview -p chr1:156000000 chr1.2mb.bam ref/hg37d5.chr1.fa
```

#### 2.4 Variants calling
A simple tool `longshot` can be used to call variants from long-reads aligned bam file. To use this tool, you need to `conda activate longshot` to activate the virtual environment.

Then,
```
longshot --bam chr1.2mb.bam --ref ref/hg37d5.chr1.fa --out mp2.longshot.vcf
``` 
will generate called variants and save in `mp2.longshot.vcf` in a vcf format. 

VCF format is plain text, and you can use `less mp2.longshot.vcf` to see what is inside this file.

#### 2.5 Statistics of called variants
We also want to see how many snps and indels are in the vcf files. To do that please run `conda deactivate` first.
```
bgzip mp2.longshot.vcf
bcftools index -f mp2.longshot.vcf.gz
bcftools stats mp2.longshot.vcf.gz | grep "TSTV"
```
will tell you 
```
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       1748    532     3.29    1748    532     3.29
```
Where TS is more than 3 times than TV.

```
bcftools stats mp2.longshot.vcf.gz | grep "^SN"
```
will tell that 
```
SN      0       number of samples:      1
SN      0       number of records:      2280
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 2280
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
Where there are 2280 snps and no indels found. 

#### 2.5.1. Evaluation of called variants.
We also validate the called variants with the gold-standard variants. 
```
bcftools isec mp2.longshot.vcf.gz vcf/vcf.chr1.2mb.vcf.gz -p longshot_perf
bcftools stats longshot_perf/0002.vcf | grep "^SN"
```
And then, we can know that 
```
SN      0       number of samples:      1
SN      0       number of records:      630
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 630
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
This overlapping suggests the recall is 630/829=0.76, but the precision is 630/2280=0.276. 


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
SN      0       number of records:      492
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 492
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```
There are much less intersected variants called from short reads and long reads. In this example, no indels and more snps (>10% of 1 millions bp region) are called, but the overlapped snps are less than 500. Maybe, more complicated variant calling tools would be used, such as `deepvariant` developed by google. One can try `deepvariant` by himself.


## After the tutorial

To do other tutorial, you might need to run `conda deactivate` to go back to the base environment for other projects if you do not practice `2.5`. ***If you still have issues to run other projects, please re-login.***

