
## Before tutorial
We have prepared some tools in the conda environment. To use the genome assembly tools in the exercise, you need to run 

```
. /shared/miniconda3/etc/profile.d/conda.sh
 conda activate base
 ```
 if you did not activate conda base.


## short read assembly
The short-read data is from 1Mb region in chr22, and the tool in this tutorial is `velvet`.

`mkdir -p ~/project/assembly/short-reads` to prepare a folder the assembly.

`cd ~/project/assembly/short-reads` to go to the assembly folder

`ln -s /shared/data/NA12878_short_30X data` to link data folder

Then we run velvet, a short read assembler using the de bruijn graph described in class.

```
velveth chr22_vel 31 -short -separate -fastq data/chr22.1mb_1.fq data/chr22.1mb_2.fq
velvetg chr22_vel
```

The results of the assembly can be found in the files 

    `chr22_vel/contigs.fa`: contig file


## long read assembly
The long-read data was sequenced for lambda, and the data is in `/shared/data/lambda_100x/lambda_100x.fastq`. In this tutorial, two tools can be used to conduct long read sequence assembly.

### Canu

First, we create a separate directory to store results.

`mkdir -p ~/project/assembly/long-reads-lambda/canu` to create an assembly folder.

Then enter the directory: `cd ~/project/assembly/long-reads-lambda/canu`


To do the assembly, run `nohup time canu -p lambda_canu -d lambda_canu genomeSize=50k -nanopore-raw /shared/data/lambda_100x/lambda_100x.fastq gnuplotTested=true useGrid=False minReadLength=500 > canu.nohup.log &`
 It would take several minutes to be done. 
 
 The results can be found in `lambda_canu/lambda_canu.contigs.fasta`

### wtdbg2

First, we create a separate directory to store results.


`mkdir -p ~/project/assembly/long-reads-lambda/wtdbg2` to create an assembly folder.
`cd ~/project/assembly/long-reads-lambda/wtdbg2`

There are several steps to use wtdbg2 for genome assembly.

Run `/shared/tools/wtdbg2/wtdbg2 -x ont -g 50k -i /shared/data/lambda_100x/lambda_100x.fastq -t 2 -fo lambda_wtdbg2` to assemble long reads. 


Run `/shared/tools/wtdbg2/wtpoa-cns -t 2 -i lambda_wtdbg2.ctg.lay.gz -fo lambda_wtdbg2.ctg.slf.fa` to derive consensus.

Additionally, the assembled sequence can be polished using the commands below
```
minimap2 -t 2 -ax map-ont lambda_wtdbg2.ctg.slf.fa reads.fa.gz | samtools sort -@2 >lambda_wtdbg2.ctg.map.srt.bam
samtools view lambda_wtdbg2.ctg.map.srt.bam | /shared/tools/wtdbg2/wtpoa-cns -t 2 -d lambda_wtdbg2.ctg.slf.fa -i - -fo lambda_wtdbg2.ctg.lrp.fa
```

Finally, the result sequence can be found in `lambda_wtdbg2.ctg.slf.fa`, and the polished sequence can be found in `lambda_wtdbg2.ctg.lrp.fa`.


### evaluation of assembly results

We next evaluate the results of the assembly. For both canu and wtdbg, only one contig is generated from the assembly of 100X data, unlike short-read assembly where hundreds of contigs are generated. The reference genome is ~48kb, so we compare the assembly with the reference genome first to evaluate completeness and accuracy.



## After tutorial

To do other tutorial, you need to run `conda deactivate` to go back to the base environment for other projects. ***If you still have issues to run other projects, please re-login.***


