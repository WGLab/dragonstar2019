mkdir -p ~/project/phenotype
cd ~/project/phenotype
disease_annotation.pl alzheimer -p -ph -logistic -out ex1
disease_annotation.pl "alzheimer;brain" -p -ph -logistic -out ex2
disease_annotation.pl /shared/tools/phenolyzer/example_phenotype.txt -f -p -ph -logistic -out ex3
disease_annotation.pl /shared/tools/phenolyzer/example_disease.txt -f -p -logistic -out ex4
cd ~/project
