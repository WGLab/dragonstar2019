## Detection of copy number variants from DNA microarrays

1. Creaet a folder for this project, using the `mkdir -p ~/project/penncnv/` command.

2. Our initial exercise will focus on examples included in the PennCNV package. 

Do `cd ~/project/penncnv` to enter the directory.

We will set up a symbolic link to make the commands below easier to write:

```
ln -s /shared/tools/PennCNV/example/example.hmm
ln -s /shared/tools/PennCNV/example/example.pfb 
ln -s /shared/tools/PennCNV/example/father.txt 
ln -s /shared/tools/PennCNV/example/mother.txt 
ln -s /shared/tools/PennCNV/example/offspring.txt 
```

Essentially, this command creates a few symbolic links to the files /shared/tools/PennCNV/example directory in the current directory. You can type a `ls -l` to see details.

3. We will first perform CNV calling on three signal intensity files to see how the procedure works:

```
detect_cnv.pl -test -hmm example.hmm -pfb example.pfb -conf -log ex1.log -out ex1.rawcnv father.txt mother.txt offspring.txt
```

This command generates CNV calls for three samples (stored in the father.txt, mother.txt and offpsring.txt files).

We can take a look of the output file to see what it contains:

```
cat ex1.rawcnv
```

As you can see, the file contains 4 CNV calls in father.txt, 2 CNV calls in mother.txt and 4 CNV calls in offspring.txt.

4. We next want to filter the CNV calls and only keep more confident calls (that contain more than 10 markers and are longer than 50kb)

```
filter_cnv.pl -numsnp 10 -length 50k -type del ex1.rawcnv > ex1.filtercnv
```

Now, take a look at the new file:

```
cat ex1.filtercnv
```

As you can see, the file contains 1 CNV call in father.txt, 1 CNV call in mother.txt and 2 CNV calls in offspring.txt. One of the CNV call is listed below:

```
chr3:3957986-4054960          numsnp=50     length=96,975      state2,cn=1 offspring.txt startsnp=rs11716390 endsnp=rs17039742 conf=219.721
```

This is a very large CNV call that is observed in offspring but not in parents. It is potentially a de novo CNV, but we want to double check this.


5. Next, we are going to infer SNP alleles and determine the de novo status of CNV calls

```
infer_snp_allele.pl -pfb example.pfb -hmm example.hmm -denovocn 1 -start rs11716390 -end rs17039742 -out ex14.geno father.txt mother.txt offspring.txt
```

This command generates a new geno file that contains the genotype calls between two SNPs, and assess whether a CNV call is likely to be a de novo CNV.

6. Finally, we will generate some plots to visually inspect the CNV calls

```
visualize_cnv.pl -format plot -signal offspring.txt ex1.filtercnv
```

You will see messages like below:

```
NOTICE: Signal values for 2 CNV regions are found in offspring.txt
NOTICE: Processing sample offspring.txt CNV chr20:10511631-10583260 with copy number of 1 ...  written to offspring.txt.chr20.10511631.JPG
NOTICE: Processing sample offspring.txt CNV chr3:3957986-4054960 with copy number of 1 ...  written to offspring.txt.chr3.3957986.JPG
```

Then run a `ls` command, and you will see two new JPG files in the directory. Try to use SFTP to transfer the two JPG files to your local computer, and examine these two image files.


## Detection of structural variants (SVs) from short-read sequencing data

We can use Delly (https://github.com/dellytools/delly) to call SVs from short-read sequencing data

1. Create a folder for this project, using the `mkdir -p ~/project/short_reads_sv/` command. Do `~/project/short_reads_sv/` to enter the directory.

2. We will use the bam files generated in the alignment training section as input. You should alreadly have the following files if you have completed the alignment section.

```
/home/biouser/project/alignment/short-reads/chr1.2mb.bwa.bam # aligner: BWA-MEM
/home/biouser/project/alignment/short-reads/chr1.2mb.mp2.bam # aligner: Minimap2
```

You can use the following commands to run Delly: 

```
delly call -g /shared/data/ref_hg37_chr1/ref/hg37d5.chr1.fa -o chr1.2mb.bwa.bam.delly.bcf /home/biouser/project/alignment/short-reads/chr1.2mb.bwa.bam
delly call -g /shared/data/ref_hg37_chr1/ref/hg37d5.chr1.fa -o chr1.2mb.mp2.bam.delly.bcf /home/biouser/project/alignment/short-reads/chr1.2mb.mp2.bam
```

In the above commands, `-g` specifies the reference fasta file and `-o` specifies the output BCF file. BCF is the binary version of VCF.

3. The program will run for a few minutes and we can get the following output files. 

```
chr1.2mb.bwa.bam.delly.bcf
chr1.2mb.bwa.bam.delly.bcf.csi
chr1.2mb.mp2.bam.delly.bcf
chr1.2mb.mp2.bam.delly.bcf.csi
```
The `.csi` file is the index file of the `.bcf` files. 

We can use bcftools to convert the `.bcf` files to VCF format: 

```
bcftools view chr1.2mb.bwa.bam.delly.vcf > chr1.2mb.bwa.bam.delly.vcf
bcftools view chr1.2mb.mp2.bam.delly.bcf > chr1.2mb.mp2.bam.delly.vcf
```

4. The know SV is 1:156526705-156528935 (deletion). We can use the following command to check if the deletion is in the vcf files: 

```
grep 156526 chr1.2mb.bwa.bam.delly.vcf chr1.2mb.mp2.bam.delly.vcf 
```

`grep` is a command to match string in text files. In this case, we want to find string `156526`. The position in the vcf file may be a few hundred base pair away from the true position so we don't grep the exact position `156526705`. 

We can get the following results:

```
chr1.2mb.bwa.bam.delly.vcf:1	156526704	DEL00000046	T	<DEL>	.	PASS	PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.7.8;CHR2=1;END=156528936;PE=101;MAPQ=60;CT=3to5;CIPOS=-21,21;CIEND=-21,21;INSLEN=0;HOMLEN=20;SR=10;SRQ=0.993827;CONSENSUS=TAAGTGTTCAGGAAGAAAAGGGGCTGGGTTGCTTTAACAAGAGGCTCTGTAAGAAGCAATTTGTCAGGCCTAGAAATTGAGTAGCTCAGCATGTAACACAGAGTGGCTGTCATGGCAGAGGGTGAGTTCCTAAGGTGGTGAGCACAAGATTGACAGGTGGCT;CE=1.94414	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:R1/1:-235.387,-19.8551,0:10000:PASS:365:141:285:0:0:102:0:66
chr1.2mb.mp2.bam.delly.vcf:1	156526704	DEL00000110	T	<DEL>	.	PASS	PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.7.8;CHR2=1;END=156528936;PE=94;MAPQ=60;CT=3to5;CIPOS=-21,21;CIEND=-21,21;INSLEN=0;HOMLEN=20;SR=10;SRQ=0.993827;CONSENSUS=TAAGTGTTCAGGAAGAAAAGGGGCTGGGTTGCTTTAACAAGAGGCTCTGTAAGAAGCAATTTGTCAGGCCTAGAAATTGAGTAGCTCAGCATGTAACACAGAGTGGCTGTCATGGCAGAGGGTGAGTTCCTAAGGTGGTGAGCACAAGATTGACAGGTGGCT;CE=1.94414	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:R1/1:-235.287,-19.8551,0:10000:PASS:365:82:280:0:0:94:0:66
```

One call was found in each vcf file.
In the first line, `1  156526704` indicates that the first breakpoint is at chr1:156526704; `CHR2=1;END=156528936` indicates that the second breakpoint is at chr1:156528936. 

In the second line, `1	156526704` indicates that the first breakpoint is at chr1:156526704; `CHR2=1;END=156528936`indicates that the second breakpoint is at chr1:156528936. 

Therefore, Delly generated the same results from `chr1.2mb.bwa.bam` and `chr1.2mb.mp2.bam`. 





