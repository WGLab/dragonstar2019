1. go to your home folder, and `mkdir -p project/annotation/` folder

2. Our initial exercise will focus on an example included in the ANNOVAR package. 

Do `cd ~/project/annotation/`, and then `mkdir ex2` and `cd ex2`

3. We will first annotate a small VCF file to see how the procedure works:

```
table_annovar.pl /shared/tools/annovar/example/ex2.vcf /shared/tools/annovar/humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput
```

The command takes the ex2.vcf file, and then a series of annotation tasks on the VCF file. Each annotation task corresponds to one protocol and one operation (such as g, r, and f). The g, r and f represent gene-based, region-based and filter-based annotation, respectively, and the different protocol names represent different databases to use for the annotation. The final output is written to myanno.hg19_multianno.txt file, as well as myanno.hg19_multianno.vcf.

we can take a look at the myanno.hg19_multianno.txt file: it is a tab-delimited file, with the first line being the header line. We can open the file in a software such as Excel to examine it in more details.

4. Next, we will try to annotate the VCF file generated from the alignment/variant exercise described previously. We want to find the refGene, cytoBand, dbNFSP scores for non-synonymous SNPs, and the allele frequency in different ethnicity groups as recorded in the gnomAD database. The command line is below:

```
table_annovar.pl ~/project/alignment/mp2.bcftools.call.bcf /shared/tools/annovar/humandb/ -buildver hg19 -out mp2 -remove -protocol refGene,cytoBand,dbnsfp35a,gnomad211_exome -operation g,r,f,f -nastring . -vcfinput
```

5. Examine the results in the output mp2.hg19_multianno.txt file.

6. Once we finish the phenolyzer exercise, we will analyze a real exome sequencing data in the anemia.vcf file. This is generated on a patient diagnosed with hemolytic anemia.





