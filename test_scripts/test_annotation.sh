mkdir -p ~/project/annotation/
cd ~/project/annotation/

mkdir ex2
cd ex2
table_annovar.pl /shared/tools/annovar/example/ex2.vcf /shared/tools/annovar/humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput -polish

cd ..
mkdir ex3
cd ex3

table_annovar.pl /shared/data/VCF/1000G_PKLR.vcf /shared/tools/annovar/humandb/ -buildver hg19 -out pklr -remove -protocol refGene,cytoBand,dbnsfp35a,gnomad211_exome -operation gx,r,f,f -nastring . -vcfinput -polish -xref /shared/tools/annovar/example/gene_xref.txt

cd ~/project

