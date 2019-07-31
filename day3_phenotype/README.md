1. go to your home folder, and `mkdir -p project/phenotype/` folder

2. Our initial exercise will focus on using the Phenolyzer software for gene prioritization:

```
disease_annotation.pl alzheimer -prediction -phenotype -logistic -out ex1
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

When the disease/phenotype contains space, and we must enclose them by quote in command line:

```
disease_annotation.pl "Amyotrophic lateral sclerosis" -prediction -phenotype -logistic -out ex2
```

When multiple phenotype terms are available, and you can separate them by ";" in command line. Alternatively, save a few phenotype terms in a file, and then supply this file in command line:

```
disease_annotation.pl /shared/tools/phenolyzer/example_phenotype.txt -file -prediction -phenotype -logistic -out ex3
```

Sometimes, if you know the disease name already and do not want to do "phenotype expansion" step, you can omit the `-ph` argument:

```
disease_annotation.pl /shared/tools/phenolyzer/example_disease.txt -file -prediction -logistic -out ex4
```

Phenolyzer can also directly process HPO terms as phenotypes. For example, one patient was seen by a genetic counselor and the list of HPO terms were recorded in a file (one term per line). We want to analyze the file and prioritize candidate genes:

```
disease_annotation.pl /shared/tools/phenolyzer/example_hpo.txt -f -p -ph -logistic -out ex5 -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
```

Examine the output files from the above commands.




3. We will next evaluate a tool called Phenomizer that uses HPO terms to predict possible diseases and genes.

The clinical presentation of a patient is described below. Please make a possible diagnosis of the disease, using the Phenomizer tool (http://compbio.charite.de/phenomizer/). The procedure will be described in class as a demo, and students will try to reproduce the results during the exercise session.

![HPO terms](case_hpo.png "Logo Title Text 1")


4. We will also explore the Phenolyzer web server, and examine some example pages. These pages are described in class.

The web page is below (depending on Internet speed, it may take a while to load the example pages)

http://phenolyzer.wglab.org/example.php


5. HPO extraction exercise

Many “case reports” will publish clinical phenotypes for patients, usually written in free texts by genetic counselors or medical geneticists. 

Example clinical notes from a paper (http://molecularcasestudies.cshlp.org/content/2/6/a001131.long): "The proband was born to a non-consanguineous couple, who had an unremarkable pregnancy history; however, at birth a large fontanel was reported. Parents and siblings were healthy and no significant family history was reported (Figure 1). The proband met all developmental milestones, except crawling, up to his first epileptic episode, which occurred at three years of age. After this episode, he lost all speech, began exhibiting autistic behavior, and also started to have frequent generalized tonic-clonic seizures. Over time, tonic, atonic, mild clonic, complex partial, myoclonic and gelastic seizures were reported in the proband. Other developmental skills, including throwing a ball, responding to his name, feeding himself with utensils and self-care skills were lost by 4-years of age. No significant conductive hearing loss, heart abnormalities or delayed bone age were found in the proband at that age. The proband was evaluated (by G.J.L.) at eleven years of age. He presented with several neurological and craniofacial abnormalities including epilepsy, ventriculomegaly, relative macrocephaly, prominent forehead, low hairline, thick eyebrows, wide-set eyes, macrodontia of upper central incisors, and full lips. Hand and foot abnormalities included clinodactyly of the fifth digit, bilateral single transverse palmar creases, brachydactyly and flat feet (Figure 2). He also had a diagnosis of cerebral folate deficiency due to the presence of folate receptor autoantibodies."

A number of tools, such as EHR-Phenolyzer, can convert free texts into HPO terms, and then use these HPO terms for gene prioritization. Some web servers are also developed for this purpose. In this bonus exercise, try to use the web server, available at https://impact2.dbmi.columbia.edu/doc2hpo/, to generate the HPO terms for this patient. Then generate prioritized gene list on this patient using Phenolyzer (please use command line to perform this analysis).










