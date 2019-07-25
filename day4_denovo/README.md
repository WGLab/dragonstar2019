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

This is a more difficult case to study, because the mutation is a de novo mutation. For your convenience, we have downloaded the VCF files for the proband, father, mother, brother 1, brother 2 and sister, to the directory `/shared/data/VCF/SCN8A`. However, we do not know in advance whether a de novo mutation, or two mutations in a recessive gene, is causal for the disease. So all the possibilities should be considered in our analysis.
