. /shared/miniconda3/etc/profile.d/conda.sh
conda activate base

mkdir -p ~/project/alignment/
cd ~/project/alignment/
mkdir short-reads
cd short-reads

ln -s /shared/data/NA12878_short_30X data
ln -s /shared/data/ref_human/chroms ref_hg38_chr22

minimap2 -ax sr ref_hg38_chr22/chr22.fa data/chr22.1mb_1.fq data/chr22.1mb_2.fq | samtools sort | samtools view -bS - > chr22.1mb.mp2.bam
samtools index chr22.1mb.mp2.bam

bwa mem ref_hg38_chr22/chr22.fa data/chr22.1mb_1.fq data/chr22.1mb_2.fq | samtools sort | samtools view -bS - > chr22.1mb.bwa.bam

samtools index chr22.1mb.bwa.bam


bcftools mpileup -f ref_hg38_chr22/chr22.fa chr22.1mb.mp2.bam | bcftools call -mv -Ob -o mp2.bcftools.call.bcf

bcftools view mp2.bcftools.call.bcf > mp2.bcftools.call.vcf

bcftools mpileup -f ref_hg38_chr22/chr22.fa chr22.1mb.bwa.bam | bcftools call -mv -Ob -o bwa.bcftools.call.bcf
bcftools view bwa.bcftools.call.bcf > bwa.bcftools.call.vcf


cd ~/project/alignment/
mkdir long-reads
cd long-reads
ln -s /shared/data/NA12878_nanopore data
ln -s /shared/data/ref_human/chroms ref_hg38_chr22
minimap2 -ax map-ont ref_hg38_chr22/chr22.fa  data/chr22.1mb.fq | samtools sort | samtools view -bS - > chr22.1mb.bam
samtools index chr22.1mb.bam

conda activate longshot
longshot --bam chr22.1mb.bam --ref ref_hg38_chr22/chr22.fa --out mp2.longshot.vcf

conda deactivate

