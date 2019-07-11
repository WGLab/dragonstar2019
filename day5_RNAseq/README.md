## 0. Introduction
This tutuorial is about how to use existing tools for RNA-seq data analysis. 
Because of the heavy tasks of aligning RNA-seq data against a reference genome. 
In this tutorial, we assume that we already have the results from the alignment and gene counts.

## 1. Preparation of directories and data files.
1.1 `mkdir -p project/RNA-seq-tutorial` and `cd project/RNA-seq-tutorial` to create tutorial directory.
1.2 `ln -s /shared/data/RNA-seq-tutorial-data data` to link the data of gene counts for further analysis.

## 2. Prepare R and data
The data analysis for differential genes will be conducted `DeSeq2` in R. So, we will open the R and load `DeSeq2` and data for next step.
Simple type `R` and `Enter` to enter R software, and then under R,
```
library("DESeq2")
counts <- read.table("NB_vs_GBM.txt",header=TRUE,stringsAsFactors=TRUE)
rownames(counts) <- counts$Geneid
counts <- counts [ ,-1]
```

To see how the data look like, 
```
head(counts)
```
will outputs
```
      Geneid BJ024 BJ028 BJ030 PC112 PC123 PC124
1    DDX11L1     0     0     0     0     2     2
2     WASH7P   114   111   165    98   161   162
3 MIR1302-10   112    58    41   151   446   397
4    FAM138A     0     0     0     0     0     0
5    OR4G11P     0     0     0     0     0     0
6      OR4F5     0     0     0     0     0     0
```

## 3. Use `DeSeq2` for data analysis
3.1 Set the group conditions
```
columndata = data.frame(row.names = colnames(counts), condition = c("NB","NB","NB","GBM","GBM","GBM"), libtype = c("single-end","single-end","single-end","single-end","single-end","single-end"))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = columndata, design =~ condition)
dds$condition <- factor(dds$condition, levels=c("NB","GBM"))
```

3.2 Get the results from `DeSeq2`
```
dds <- DESeq(dds)
res <- results(dds)
sortres <- res[order(res$padj),]
summary(res)
```
It might take several seconds, and then, the last line will output the summary of the results
```
out of 33296 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 4699, 14%
LFC < 0 (down)   : 4767, 14%
outliers [1]     : 191, 0.57%
low counts [2]   : 10049, 30%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

## 4. Plot the results and save them
4.1 Plot the results.
```
pdf("m_plotMA.pdf");
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()
```
The plot will be saved in `m_plotMA.pdf`.

4.2 Plot gene read counts
```
pdf("m_plotCounts.pdf");
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()
```

4.3 Plot PCA
```
rld <- rlog(dds, blind=FALSE)
pdf("m_plotPCA.pdf");
plotPCA(rld, intgroup=c("condition", "libtype"))
dev.off()
```

4.4 Using `ggplot`
```
data <- plotPCA(rld, intgroup=c("condition", "libtype"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pdf("m_plotPCA_ggplot.pdf");
ggplot(data, aes(PC1, PC2, color=condition, shape=libtype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()
```

4.5 Plot p-value distribution
```
pdf("m_hist.pdf");
hist(res$pvalue, br=20)
dev.off()
```

## 5. Save the results
```
write.csv(as.data.frame(sortres),file="NB_v_GBM.csv")
```
And you will find `NB_v_GBM.csv`.

After that, one needs to exit R by `Ctrl + D` and then `n` for `Save workspace image? [y/n/c]:`. And then,
```
awk -F ',' '{if(($3 > 1 || $3 < -1) && length($1)>2) print $1}' NB_v_GBM.csv  | head -101 | sed 's/\"//g' > NB_v_GBM.csv.top100_genes.txt
```
The top 100 genes with folder change > 2 (or < 0.5) will be output to `NB_v_GBM.csv.top100_genes.txt`. One can copy the gene list and paste to Enrichr for Gene Ontology, pathway, and TF-target enrichment analyses.


