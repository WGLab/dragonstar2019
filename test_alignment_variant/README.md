To use installed softwares, you need to run `. /shared/miniconda3/etc/profile.d/conda.sh` and then `conda activate base`.

1. go to your home folder, and `mkdir -p project/alignment/` folder

2. Alignment and variant calling based on short reads
* 2.1 `cd ~/project/alignment/`, and then `mkdir short-reads` and `cd short-reads`
* 2.2 Link data by `ln -s /shared/data/NA12878_short_30X data`
* 2.3 Link reference genome by `ln -s /shared/data/ref_human/chroms ref_hg38_chr22`
* 2.4 Map short reads to a reference genome. In this example, the reads are from a 1MB region in chr22. There are two ways to do the alignment.

    * 2.4.1 Alignment with minimap2. 
`minimap2 -ax sr ref_hg38_chr22/chr22.fa data/chr22.1mb_1.fq data/chr22.1mb_2.fq | samtools sort | samtools view -bS - > chr22.1mb.mp2.bam`
Reads will be aligned with chr22, and then sorted and saved into a bam file.
After that, `samtools index chr22.1mb.mp2.bam` is used to build index, and `*.bai` will created for the bam file.

    * 2.4.2 Alignment with bwa-mem. 
`bwa mem ref_hg38_chr22/chr22.fa data/chr22.1mb_1.fq data/chr22.1mb_2.fq | samtools sort | samtools view -bS - > chr22.1mb.bwa.bam`. And then index the bam using `samtools index chr22.1mb.bwa.bam`

* 2.5 One can see what is bam in two ways.
    * 2.5.1 `samtools view chr22.1mb.bwa.bam | less` to see how to represent alignment for each reads. `type q` to exit.
    * 2.5.2 `samtools tview -p chr22:25499651 chr22.1mb.bwa.bam ref_hg38_chr22/chr22.fa` to see how alignments parallel with the reference genome.

* 2.6 A simple way for variant calling is to use bcftools. We can generate variant calling for both bam files generatd above.
    * 2.6.1 for `chr22.1mb.mp2.bam`. 
`bcftools mpileup -f ref_hg38_chr22/chr22.fa chr22.1mb.mp2.bam | bcftools call -mv -Ob -o mp2.bcftools.call.bcf` will call variants and save to `mp2.bcftools.call.bcf` in a bcf format. You can view the file using `bcftools view mp2.bcftools.call.bcf | less` or `bcftools view -i '%QUAL>=20' mp2.bcftools.call.bcf | less` to only view those higher quality (>20) variants.
`bcftools view mp2.bcftools.call.bcf > mp2.bcftools.call.vcf` can convert bcf to vcf;

    * 2.6.2 for `chr22.1mb.bwa.bam`
Simiarly, `bcftools mpileup -f ref_hg38_chr22/chr22.fa chr22.1mb.bwa.bam | bcftools call -mv -Ob -o bwa.bcftools.call.bcf` to call variants and save to `bwa.bcftools.call.bcf`. You can also use `bcf view` to what variants have been called.
`bcftools view bwa.bcftools.call.bcf > bwa.bcftools.call.vcf` can convert bcf to vcf for further analysis.

3. Alignment and variant calling based on long reads
    * 3.1 `cd ~/project/alignment/`, and then `mkdir long-reads` and `cd long-reads`
    * 3.2 Link data by `ln -s /shared/data/NA12878_nanopore data`
    * 3.3 Link reference genome by `ln -s /shared/data/ref_human/chroms ref_hg38_chr22`
    * 3.4 Map short reads to a reference genome. 
      ```
      minimap2 -ax map-ont ref_hg38_chr22/chr22.fa  data/chr22.1mb.fq | samtools sort | samtools view -bS - > chr22.1mb.bam
      samtools index chr22.1mb.bam
      ```
      The commands above will align long reads in `data/chr22.1mb.fq`, and then sort/save alignment into a bam file `index chr22.1mb.bam`. An index is also created so that you can use `samtools view` or `samtools tview` to view the alignment in the bam file.

    * 3.5 A simple tool `longshot` can be used to call variants from long-reads aligned bam file. To use this tool, you need to `conda activate longshot` to activate the virtual environment.
   Then, `longshot --bam chr22.1mb.bam --ref ref_hg38_chr22/chr22.fa --out mp2.longshot.vcf` will generate called variants and save in `mp2.longshot.vcf` in a vcf format. 

   VCF format is plain text, and you can use `less mp2.longshot.vcf` to see what is inside this file.

To do other tutorial, you might need to run `conda deactivate` to go back to the base environment.

