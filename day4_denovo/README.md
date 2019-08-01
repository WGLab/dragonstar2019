# Identification of candidate genes for a patient with hemolytic anemia

In this exercise, we will analyze an exome file in VCF format for a patient and find the likely causal gene and a possible genetic diagnosis.

The VCF file is at /shared/data/vcf/anemia.vcf. The phenotype is 'hemolytic anemia'.

We will perform variant annotation first and generate multianno file.

The student can write a script themselves to filter the multianno file, or use Excel to filter the file. (the `variants_reduction.pl` program in ANNOVAR can also be used.) The goal is to reduce the search space to only a handful of likely disease-causal variants or disease-causal genes, to faciliate genetic diagnosis.

At the same time, run phenolyzer to prioritze genes (hint: use "hemolytic anemia" as the phenotype).

Then combine the genotype/phenotype results together to find disease causal gene (PKLR).

This is an analysis task that the students should try to solve the exome case themselves without following commands. Therefore, we do not provide example command here.

# Identification of causal de novo mutation for a child presenting with seizures and developmental delays

This is a case that we published a few years ago, on finding disease genes on a child with seizures and developmental delays. In the original publication, we provided detailed genotype information and phenotype information (including movies and photos and clinical phenotype descriptions) on the patient. These files are available in the [supplementary materials](http://molecularcasestudies.cshlp.org/content/early/2016/07/19/mcs.a001073/suppl/DC1).

The original publication used several web servers, including the Phenolyzer server and the wANNOVAR server, to identify the disease causal gene. Your challenge here is not to use any external web servers! Instead, try to use only command line tools to solve this case. You may want to quickly browser through the [paper](http://molecularcasestudies.cshlp.org/content/2/6/a001073.full) first to get an idea of the family and the phenotypic presentation.

This is a more difficult case to study, because the mutation is a de novo mutation. For your convenience, we have downloaded the VCF files for the proband, father, mother, brother 1, brother 2 and sister, to the directory `/shared/data/VCF/SCN8A`. 

To make things easier, below is a list of clinical phenotypes on the patient as reported in the paper. You can save it in a file.

```
HP:0200134
HP:0010818
HP:0007193
HP:0002353
HP:0001263
HP:0006834
HP:0002376
HP:0001344
HP:0010864
HP:0001270
HP:0001290
HP:0012389
HP:0009062
HP:0000467
HP:0002063
HP:0001257
HP:0001531
HP:0011471
HP:0002020
HP:0002015
HP:0002880
HP:0012418
HP:0100765
HP:0002870
HP:0000248
HP:0000337
HP:0000431
HP:0000430
HP:0000293
HP:0000212
HP:0000347
HP:0002267
HP:0002345
HP:0000643
HP:0001283
HP:0000639
HP:0001347
HP:0008763
```

You can first do the analysis by looking at predicted de novo mutations in the proband (i.e. those mutations in proband that are not in parents).

Hint: since the clinical phenotypes are well characterized already, a command such as `disease_annotation.pl pheno_list.txt -file -logistic -out hpo2gene -phenotype -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25` can generate candidate genes for these well characterized phenotypes. (If phenotypes may have error or incomplete, then you should add -prediction in the command line to expand the phenotype to similar phenotype terms)

Hint: since we want to find predicted de novo mutations in the proband, we can use `bcftools isec` on all the VCF files to find mutations only in the first file but not other family members, through the `-C` argument. Then do annotation on this file, and use `bcftools view` to filter and generate the subset of predicted de novo mutations that are likely pathogenic (non-synonymous, nonsense, frameshift, splicing, and not observed in gnomad database). An expression like `(ExonicFunc.refGene="nonsynonymous_SNV" || ExonicFunc.refGene="frameshift_deletion" || ExonicFunc.refGene="frameshift_insertion" || ExonicFunc.refGene="frameshift_block_substitution" || ExonicFunc.refGene="stopgain" || ExonicFunc.refGene="stoploss" || Func.refGene="splicing") && (AF_popmax="." || AF_popmax<0.0001) && FORMAT/DP>=20` may be used.


Since we do not know in advance whether a de novo mutation, or two deleterious mutations (one inherited from father, one inherited from mother) in a recessive gene, is causal for the disease. So next we should consider this possibility.




