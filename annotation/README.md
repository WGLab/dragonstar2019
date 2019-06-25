1. go to your home folder, and `mkdir -p ~/project/annotation/` folder

2. Our initial exercise will focus on an example included in the ANNOVAR package. 

Do `cd ~/project/annotation/`, and then `mkdir ex2` and `cd ex2`

3. We will first annotate a small VCF file to see how the procedure works:

```
table_annovar.pl /shared/tools/annovar/example/ex2.vcf /shared/tools/annovar/humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp35a -operation g,r,f,f,f -nastring . -vcfinput -polish
```

The command takes the `ex2.vcf` file, and then a series of annotation tasks on the VCF file. Each annotation task corresponds to one protocol and one operation (such as `g`, `r`, and `f`). The `g`, `r` and `f` represent gene-based, region-based and filter-based annotation, respectively, and the different protocol names represent different databases to use for the annotation. The final output is written to `myanno.hg19_multianno.txt` file, as well as `myanno.hg19_multianno.vcf`.

we can take a look at the `myanno.hg19_multianno.txt` file: it is a tab-delimited file, with the first line being the header line. We can open the file in a software such as Excel to examine it in more details.

4. Next, we will try to annotate the VCF file generated from the exercise described previously. We want to find the refGene, cytoBand, dbNFSP scores for non-synonymous SNPs, and the allele frequency in different ethnicity groups as recorded in the gnomAD database. We will only use the exome subset of the gnomAD data, since the genome subset of gnomAD is too large to be used in this exercise.

We first go up one directory:

```
cd ..
mkdir ex3
cd ex3
```

Then perform the annotation. The command line is below:

```
table_annovar.pl /shared/data/VCF/1000G_PKLR.vcf /shared/tools/annovar/humandb/ -buildver hg19 -out pklr -remove -protocol refGene,cytoBand,dbnsfp35a,gnomad211_exome -operation gx,r,f,f -nastring . -vcfinput -polish -xref /shared/tools/annovar/example/gene_xref.txt
```

The `-operation` argument tells ANNOVAR which operations to use for each of the protocols: `g` means gene-based, `gx` means gene-based with cross-reference annotation (from `-xref` argument), `r` means region-based and `f` means filter-based. 

In the command above, we used `-xreffile` argument to attach additional annotation to gene names. If the file contains header line, it is possible to provide mulitple pieces of annotations to genes (rather than just one single column). To illustrate this, we can check the first 2 lines (including the header line) of the `/shared/tools/annovar/example/gene_fullxref.txt` file:

```
head -n 2 /shared/tools/annovar/example/gene_fullxref.txt
```

The results are shown below:
```
#Gene_name      pLi     pRec    pNull   Gene_full_name  Function_description    Disease_description     Tissue_specificity(Uniprot)     Expression(egenetics)  Expression(GNF/Atlas)    P(HI)   P(rec)  RVIS    RVIS_percentile GDI     GDI-Phred
A1BG    9.0649236354772e-05     0.786086131023045       0.2138232197406 alpha-1-B glycoprotein  .       .       TISSUE SPECIFICITY: Plasma.;    unclassifiable (Anatomical System);amygdala;prostate;lung;islets of Langerhans;liver;spleen;germinal center;brain;thymus;       fetal liver;liver;fetal lung;trigeminal ganglion;       0.07384 0.31615 -0.466531444    23.51380042     79.3774 1.88274
```

The header line starts with #. The cross-reference file then contains 15 types of annotations for genes.

5. Examine the results in the output `pklr.hg19_multianno.txt` file.

6. Once we finish the `test_phenotype` exercise, we will analyze a real exome sequencing data in the `anemia.vcf` file (this is generated on a patient diagnosed with hemolytic anemia). See the `test_exome` for more details on the exercise.





