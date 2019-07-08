
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

