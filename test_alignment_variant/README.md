## Before the tutorial

1. To use installed softwares, you need to run `. /shared/miniconda3/etc/profile.d/conda.sh` and then `conda activate base`.
2. go to your home folder, and `mkdir -p project/alignment/` folder to create an alignment folder

## The tutorial for read alignment and variants calling

### 1. Short read alignment and variants calling
#### 1.1 Preparation of the folder and data
1. `cd ~/project/alignment/`, and then `mkdir short-reads` and `cd short-reads`
2. Link data by `ln -s /shared/data/NA12878_short_30X data`
3. Link reference genome by `ln -s /shared/data/ref_human/chroms ref_hg38_chr22`

#### 1.2 Short reads alignment
In this tutorial, the reads are from a 1MB region in chr22. There are two ways to do the alignment.

##### 1.2.1. Alignment with minimap2. 
```
minimap2 -ax sr ref_hg38_chr22/chr22.fa data/chr22.1mb_1.fq data/chr22.1mb_2.fq | samtools sort | samtools view -bS - > chr22.1mb.mp2.bam
```
Reads will be aligned with chr22, and then sorted and saved into a bam file.
After that, `samtools index chr22.1mb.mp2.bam` is used to build index, and `*.bai` will created for the bam file.

##### 1.2.2. Alignment with bwa-mem. 
```
bwa mem ref_hg38_chr22/chr22.fa data/chr22.1mb_1.fq data/chr22.1mb_2.fq | samtools sort | samtools view -bS - > chr22.1mb.bwa.bam
```
And then index the bam using `samtools index chr22.1mb.bwa.bam`

#### 1.3 View bam files
One can see what is bam in two ways.

##### 1.3.1 View bam records
```
samtools view chr22.1mb.bwa.bam | less
``` 
to see how to represent alignment for each reads. 
`type q` to exit.

###### 1.3.2 tview bam files
```
samtools tview -p chr22:25499651 chr22.1mb.bwa.bam ref_hg38_chr22/chr22.fa
``` 
to see how alignments parallel with the reference genome.

#### 1.4 Variants calling
A simple way for variant calling is to use bcftools. We can generate variant calling for both bam files generatd above.

1. For `chr22.1mb.mp2.bam`. 
```
bcftools mpileup -f ref_hg38_chr22/chr22.fa chr22.1mb.mp2.bam | bcftools call -mv -Ob -o mp2.bcftools.call.bcf
```
will call variants and save to `mp2.bcftools.call.bcf` in a bcf format. 
You can view the file using 
```
bcftools view mp2.bcftools.call.bcf | less
``` 
or 
```bcftools view -i '%QUAL>=20' mp2.bcftools.call.bcf | less
``` 
to only view those higher quality (>20) variants.

If you want to need vcf files for analysis, please use
```
bcftools view mp2.bcftools.call.bcf > mp2.bcftools.call.vcf
``` 
to convert bcf to vcf;

2. For `chr22.1mb.bwa.bam`
Simiarly, 
```
bcftools mpileup -f ref_hg38_chr22/chr22.fa chr22.1mb.bwa.bam | bcftools call -mv -Ob -o bwa.bcftools.call.bcf
``` 
to call variants and save to `bwa.bcftools.call.bcf`. 
You can also use `bcf view` to what variants have been called.

Also, 
```
bcftools view bwa.bcftools.call.bcf > bwa.bcftools.call.vcf
``` 
can be used to convert bcf to vcf for further analysis.

### 2. Long read alignment and variants calling
#### 2.1 Preparation of the folder and data
1. `cd ~/project/alignment/`, and then `mkdir long-reads` and `cd long-reads`
2. Link data by `ln -s /shared/data/NA12878_nanopore data`
3. Link reference genome by `ln -s /shared/data/ref_human/chroms ref_hg38_chr22`

#### 2.2 Long reads alignment
```
minimap2 -ax map-ont ref_hg38_chr22/chr22.fa  data/chr22.1mb.fq | samtools sort | samtools view -bS - > chr22.1mb.bam
samtools index chr22.1mb.bam
```
The commands above will align long reads in `data/chr22.1mb.fq`, and then sort/save alignment into a bam file `index chr22.1mb.bam`. 
An index is also created so that you can use `samtools view` or `samtools tview` to view the alignment in the bam file.

#### 2.3 View bam files
One can try one of the commands below to view the bam files.
```
samtools view chr22.1mb.bam | less
samtools tview -p chr22:25499651 chr22.1mb.bam ref_hg38_chr22/chr22.fa
```

#### 2.4 Variants calling
A simple tool `longshot` can be used to call variants from long-reads aligned bam file. To use this tool, you need to `conda activate longshot` to activate the virtual environment.

Then,
```
longshot --bam chr22.1mb.bam --ref ref_hg38_chr22/chr22.fa --out mp2.longshot.vcf
``` 
will generate called variants and save in `mp2.longshot.vcf` in a vcf format. 

VCF format is plain text, and you can use `less mp2.longshot.vcf` to see what is inside this file.

## After the tutorial

To do other tutorial, you might need to run `conda deactivate` to go back to the base environment for other projects.

