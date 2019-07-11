
## Before tutorial
We have prepared some tools in the conda environment. To use the genome assembly tools in the exercise, you need to run 

```
. /shared/miniconda3/etc/profile.d/conda.sh
 conda activate base
 ```
 if you did not activate conda base.


## 1. Short read assembly
The short-read data is from simulated short-reads for lambda, and the tool in this tutorial is `velvet`.

### 1.1 Preparation of folders and data
1. `mkdir -p ~/project/assembly/short-reads` to prepare a folder the assembly.
2. `cd ~/project/assembly/short-reads` to go to the assembly folder
3. `ln -s /shared/data/lambda_short data` to link data folder

### 1.2 Running `velvet`
Then we run `velvet`, a short read assembler using the de bruijn graph described in class.

```
velveth lambda_vel 31 -short -separate -fastq data/art.lambda.cov1001.fq data/art.lambda.cov1002.fq
velvetg lambda_vel
```

#### 1.3 Assembly results
The results of the assembly can be found in the files 

`lambda_vel/contigs.fa`: contig file

In this file, one can see 17 contig in this file

## 2. Long read assembly
The long-read data was sequenced for lambda, and the data is in `/shared/data/lambda_100x/lambda_100x.fastq`. In this tutorial, two tools can be used to conduct long read sequence assembly.

### 2.1 Canu
#### 2.1.1 Preparation of folders and data
First, we create an assembly folder to store results.
```
mkdir -p ~/project/assembly/long-reads-lambda/canu
```
Then enter the directory: `cd ~/project/assembly/long-reads-lambda/canu`

#### 2.1.2 Running `canu`
To do the assembly, run 
```
nohup time canu -p lambda_canu -d lambda_canu genomeSize=50k -nanopore-raw /shared/data/lambda_100x/lambda_100x.fastq gnuplotTested=true useGrid=False minReadLength=500 > canu.nohup.log &
```
It would take several minutes to be done. 
 
#### 2.1.3 Assembly results
The results can be found in `lambda_canu/lambda_canu.contigs.fasta`

### 2.2 wtdbg2
#### 2.2.1 Preparation of folders and data
First, we create an assembly folderto store results.
```
mkdir -p ~/project/assembly/long-reads-lambda/wtdbg2
cd ~/project/assembly/long-reads-lambda/wtdbg2
```

#### 2.2.2 Run `wtdbg2` step-by-step
There are several steps to use `wtdbg2` for genome assembly.

```
/shared/tools/wtdbg2/wtdbg2 -x ont -g 50k -i /shared/data/lambda_100x/lambda_100x.fastq -t 2 -fo lambda_wtdbg2
```
to assemble long reads. 

```
/shared/tools/wtdbg2/wtpoa-cns -t 2 -i lambda_wtdbg2.ctg.lay.gz -fo lambda_wtdbg2.ctg.slf.fa
```
to derive consensus.

#### 2.2.3 Polishing
Additionally, the assembled sequence can be polished using the commands below
```
minimap2 -t 2 -ax map-ont lambda_wtdbg2.ctg.slf.fa /shared/data/lambda_100x/lambda_100x.fastq | samtools sort -@2 >lambda_wtdbg2.ctg.map.srt.bam
samtools view lambda_wtdbg2.ctg.map.srt.bam | /shared/tools/wtdbg2/wtpoa-cns -t 2 -d lambda_wtdbg2.ctg.slf.fa -i - -fo lambda_wtdbg2.ctg.lrp.fa
```

#### 2.2.4 Assembly results
Finally, the result sequence can be found in `lambda_wtdbg2.ctg.slf.fa`, and the polished sequence can be found in `lambda_wtdbg2.ctg.lrp.fa`.


## 3. Evaluation of assembly results

We next evaluate the results of the assembly. For both canu and wtdbg, only one contig is generated from the assembly of 100X data, unlike short-read assembly where hundreds of contigs are generated. The reference genome is ~48kb, so we compare the assembly with the reference genome first to evaluate completeness and accuracy.



## After tutorial

To do other tutorial, you need to run `conda deactivate` to go back to the base environment for other projects. ***If you still have issues to run other projects, please re-login.***


