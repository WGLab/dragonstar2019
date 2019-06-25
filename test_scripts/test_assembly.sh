. /shared/miniconda3/etc/profile.d/conda.sh
conda activate base

mkdir -p ~/project/assembly/short-reads
cd ~/project/assembly/short-reads
ln -s /shared/data/NA12878_short_30X data
velveth chr22_vel 31 -short -separate -fastq data/chr22.1mb_1.fq data/chr22.1mb_2.fq
velvetg chr22_vel

mkdir -p ~/project/assembly/long-reads-lambda/canu
cd ~/project/assembly/long-reads-lambda/canu
canu -p lambda_canu -d lambda_canu genomeSize=50k -nanopore-raw /shared/data/lambda_100x/lambda_100x.fastq gnuplotTested=true useGrid=False minReadLength=500 > canu.nohup.log

mkdir -p ~/project/assembly/long-reads-lambda/wtdbg2
cd ~/project/assembly/long-reads-lambda/wtdbg2
/shared/tools/wtdbg2/wtdbg2 -x ont -g 50k -i /shared/data/lambda_100x/lambda_100x.fastq -t 2 -fo lambda_wtdbg2
/shared/tools/wtdbg2/wtpoa-cns -t 2 -i lambda_wtdbg2.ctg.lay.gz -fo lambda_wtdbg2.ctg.slf.fa

minimap2 -t 2 -ax map-ont lambda_wtdbg2.ctg.slf.fa reads.fa.gz | samtools sort -@2 >lambda_wtdbg2.ctg.map.srt.bam
samtools view lambda_wtdbg2.ctg.map.srt.bam | /shared/tools/wtdbg2/wtpoa-cns -t 2 -d lambda_wtdbg2.ctg.slf.fa -i - -fo lambda_wtdbg2.ctg.lrp.fa

conda deactivate
