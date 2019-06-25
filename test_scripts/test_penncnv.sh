mkdir -p ~/project/penncnv

cd ~/project/penncnv


ln -s /shared/tools/PennCNV/example/example.hmm
ln -s /shared/tools/PennCNV/example/example.pfb 
ln -s /shared/tools/PennCNV/example/father.txt 
ln -s /shared/tools/PennCNV/example/mother.txt 
ln -s /shared/tools/PennCNV/example/offspring.txt 

detect_cnv.pl -test -hmm example.hmm -pfb example.pfb -conf -log ex1.log -out ex1.rawcnv father.txt mother.txt offspring.txt
filter_cnv.pl -numsnp 10 -length 50k -type del ex1.rawcnv > ex1.filtercnv
infer_snp_allele.pl -pfb example.pfb -hmm example.hmm -denovocn 1 -start rs11716390 -end rs17039742 -out ex14.geno father.txt mother.txt offspring.txt
visualize_cnv.pl -format plot -signal offspring.txt ex1.filtercnv

cd ~/project

