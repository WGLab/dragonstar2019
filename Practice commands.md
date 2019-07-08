### Preparation

1. To use installed softwares, you need to run 
```
. "/shared/miniconda3/etc/profile.d/conda.sh"
conda activate base
```

You will see that your command line prompt is now changed to something that look like `(base) [biouser@test-0010 ]$ `. This indicates that you are already in the base environment in conda, and we can proceed with the exercise.

2. Prepare folders
2.1 Go to your home folder (you can type `cd` to automatically go to your home folder, and `mkdir -p ~/project/startup/` folder, and then `cd ~/project/startup/`. All our tests will be done in this directory.
2.2 `ln -s /shared/data/practice_data data` to link the fq, bam and vcf files for practice.

### Practice basic commands for fq and fa files

How many reads in a fq.
```
sed -n '1~4 p' data/sr.chr1.2mb_1.fq | wc
```

Get sequences for all reads
```
sed -n '1~4 p' data/sr.chr1.2mb_1.fq | less
```
or with line number
```
sed -n '2~4 p' data/sr.chr1.2mb_1.fq | awk '{print NR" : "$0}' | less
```

Get information after '+' for each reads
```
cat data/sr.chr1.2mb_1.fq | paste - - - - | awk '{print $3"\t"$4}' | less
```


Check the number of reads whose sequence containing 'N'
```
sed -n '2~4 p' data/sr.chr1.2mb_1.fq | grep "N" | wc
```

Check the total number of bases in a fq.
```
sed -n '2~4 p' data/sr.chr1.2mb_1.fq | wc
```
The last number if the total number of bases.


Count the total number of 'N' in a fq.
```
cat data/sr.chr1.2mb_1.fq | paste - - - - | cut -f 2 | grep -o [N] | wc
```

Count the number of 'ATCGNatcg' at the first 100 reads.
```
cat data/sr.chr1.2mb_1.fq  | head -n 100 | paste - - - - | cut -f 2 | grep -o [ATCGNatcg] | sort | uniq -c
```

Convert fq to fa
```
cat data/sr.chr1.2mb_1.fq | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > sr.chr1.2mb_1.fa
```

Delete 'N' from sequences of reads
```
sed -e '/^[^>]/s/N//g' sr.chr1.2mb_1.fa | less
```

Only output reads whose sequence length larger than 65.
```
cat data/sr.chr1.2mb_1.fq | paste - - | awk '{if (length($2)>65) print $1"\n"$2}' | less
```


### Practice basic commands for bam/sam files

Get basic statistics from a bam
```
samtools stats data/chr1.2mb.mp2.bam | grep ^SN | cut -f 2-
```

Check reads which failed to align in a bam
```
samtools view data/chr1.2mb.mp2.bam | awk '{if (and($2,4)) print NR" : "$0}' | less
```

Find FLAG distribution in a bam
```
samtools view data/chr1.2mb.mp2.bam | cut -f 2 | sort | uniq -c
```

Check those alignment whose mapping quality > 30.
```
samtools view data/chr1.2mb.mp2.bam | awk '{if ($5>30) print $0}' | less
```

Check the read alignment distribution according to chromosomes.
```
samtools view data/chr1.2mb.mp2.bam | cut -f 3 | sort | uniq -c
```

Find exact match in a bam
```
samtools view data/chr1.2mb.mp2.bam | cut -f 6| awk '!/[NDISH]/{print $0}' | less
```

Find alignment with insertion and deletion
```
samtools view data/chr1.2mb.mp2.bam | cut -f 6| awk '/[DI]/{print $0}' | less
```

### Practice basic commands for bcf/vcf files

Get only snps or indels
```
bcftools view -v snps data/mp2.bcftools.call.bcf | less
bcftools view -v indels data/mp2.bcftools.call.bcf | less
```

Get insertion only
```
bcftools view -v indels data/mp2.bcftools.call.bcf | grep "^[^#]" | awk '{if (length($4)<length($5)) print $0}' | less
```
and gent deletions only
```
bcftools view -v indels data/mp2.bcftools.call.bcf | grep "^[^#]" | awk '{if (length($4)>length($5)) print $0}' | less
```

Get the distribution of snps possibliity
```
bcftools view -v snps data/mp2.bcftools.call.bcf | grep "^[^#]" | awk '{print $4"->"$5}' | sort | uniq -c
```

Get averaged coverage for snps
```
bcftools view -v snps data/mp2.bcftools.call.bcf | grep "^[^#]" | cut -d "=" -f 2| cut -d ";" -f 1 | awk '{total+=$1} END {print total/NR}' | less
```

Get averaged coverage for indels
```
bcftools view -v indels data/mp2.bcftools.call.bcf | grep "^[^#]" | cut -d "=" -f 2| cut -d ";" -f 1 | awk '{total+=$1} END {print total/NR}' | less
```

Get snps and indels whose coverage > 30
```
bcftools view -i 'MIN(DP)>30' data/mp2.bcftools.call.bcf | less
```

Get snps and indels whose coverage > 30 and quality > 10
```
bcftools view -i 'QUAL>10 & MIN(DP)>30' data/mp2.bcftools.call.bcf | less
```

Get genotype not '1/1'
```
bcftools view -i 'GT!="1/1"' data/mp2.bcftools.call.bcf | less
```

## Advanced practice of processing VCF files

You can use tabix to extract subsets of the vcf files from the 1000genomes websites. Thanks to the fact that tabix uses a index file, you will be able to download only portions of the files, without having to download everything in local.

```
tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 1:155259084-155271225 > 1000G_PKLR.vcf
```

This command takes the genomic region "1:155259084-155271225" from the 1000 Genomes Project website and save it as the 1000G_PKLR.vcf file.

The input file is hosted in the FTP server contains the 1000 Genomes Project phase3 release of variant calls. This variant set contains 2504 individuals from 26 populations.

Sometimes, due to network problems, the program above does not return a result within a few minutes. In that case, do not worry about it, just stop the program (pressing "Ctrl+C" will stop the program). We already saved a copy of the output file in `/shared/data/VCF/1000G_PKLR.vcf`, and we can just use this file in the following steps.

Next, we want to extract a few samples from the VCF file.
```
bcftools view -s HG00376,NA11933,NA12282 /shared/data/VCF/1000G_PKLR.vcf > test1.vcf
```

In the output VCF file, the original INFO field from the original VCF file is still present. We want to remove the INFO field:

```
bcftools annotate -x INFO test1.vcf > test2.vcf
```

Now examine the new vcf file, you can see that the INFO field is no longer there.


```
bcftools view --min-ac=1 test2.vcf > test3.vcf
```

Now examine the new vcf file, you can see that the INFO field is updated with the AC and AN.

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00376 NA11933 NA12282
1       155260382       rs61755431      C       T       100     PASS    AC=3;AN=6       GT      1|0     1|0     0|1
1       155260780       rs576493870     AAAT    A       100     PASS    AC=3;AN=6       GT      0|1     0|1     0|1
1       155265177       rs2071053       A       G       100     PASS    AC=1;AN=6       GT      0|0     0|1     0|0
1       155266335       rs8177968       C       T       100     PASS    AC=1;AN=6       GT      0|0     0|1     0|0
1       155268120       rs12067675      T       C       100     PASS    AC=1;AN=6       GT      0|0     0|1     0|0
1       155268425       rs12741350      C       T       100     PASS    AC=1;AN=6       GT      0|0     0|1     0|0
1       155268958       rs11264357      C       T       100     PASS    AC=1;AN=6       GT      0|0     0|1     0|0
1       155269776       rs3020781       A       G       100     PASS    AC=1;AN=6       GT      0|0     0|1     0|0
```

So now we have 8 sites that are polymorphic in the VCF file (in other words, we found eight SNPs that have mutations in at least one of the three subjects.

You can make a statistics on the test3.vcf file:

```
bcftools stats test3.vcf
```

We will see the information below:

```
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      3
SN      0       number of records:      8
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 7
SN      0       number of MNPs: 0
SN      0       number of indels:       1
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
```


### For other practice
To do other tutorial, you might need to run `conda deactivate` to go back to the base environment. You might have errors if you do not deactivate the environment.

