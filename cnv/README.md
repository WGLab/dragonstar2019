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





