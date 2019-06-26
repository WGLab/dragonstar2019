mkdir -p ~/project/phenotype
cd ~/project/phenotype
disease_annotation.pl alzheimer -prediction -phenotype -logistic -out ex1
disease_annotation.pl "Amyotrophic lateral sclerosis" -prediction -phenotype -logistic -out ex2
disease_annotation.pl /shared/tools/phenolyzer/example_phenotype.txt -file -prediction -phenotype -logistic -out ex3
disease_annotation.pl /shared/tools/phenolyzer/example_disease.txt -file -prediction -logistic -out ex4
disease_annotation.pl /shared/tools/phenolyzer/example_hpo.txt -file -prediction -phenotype -logistic -out ex5
cd ~/project
