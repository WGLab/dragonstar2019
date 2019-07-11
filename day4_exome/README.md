In this exercise, we will analyze an exome file in VCF format for a patient and find the likely causal gene and a possible genetic diagnosis.

The VCF file is at /shared/data/vcf/anemia.vcf. The phenotype is 'hemolytic anemia'.

We will perform variant annotation first (see examples in test_annovar) and generate multianno file.

The student can write a script themselves to filter the multianno file, or use Excel to filter the file. (the variants_reduction.pl file can also be used)

At the same time, run phenolyzer to prioritze genes (see examples in test_phenotype).

Then combine the genotype/phenotype results together to find disease causal gene (PKLR).

This is an advanced analysis task that the students should try to solve the exome case themselves without following commands.
