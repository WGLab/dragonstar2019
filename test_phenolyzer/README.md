1. go to your home folder, and `mkdir -p project/phenotype/` folder

2. Our initial exercise will focus on using the Phenolyzer software for gene prioritization:

```
disease_annotation.pl alzheimer -p -ph -logistic -out ex1
```

The input in this case is `alzheimer`, and the output are written to a few files. The file that we are interested in is `ex1.final_gene_list`. We can check the first 10 lines of the file:

```
[biouser@ecs-3e58 test_phenolyzer]$ head ex1.final_gene_list
Rank    Gene    ID      Score   Status
1       PSEN2   5664    1       SeedGene
2       APP     351     0.9799  SeedGene
3       PSEN1   5663    0.956   SeedGene
4       APOE    348     0.414   SeedGene
5       GATA1   2623    0.3217  SeedGene
6       PRNP    5621    0.3088  SeedGene
7       CACNA1G 8913    0.3053  SeedGene
8       CASP3   836     0.2479  Predicted
9       SORL1   6653    0.2223  SeedGene
```

As you can see, a few genes are prioritized as the disease candidate genes, each with a Phenolyzer score.

3. We will next evaluate a tool called Phenomizer that uses HPO terms to predict possible diseases and genes.

The clinical presentation of a patient is described below. Please make a possible diagnosis of the disease, using the Phenomizer tool (http://compbio.charite.de/phenomizer/). 








