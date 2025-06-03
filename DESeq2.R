#Author: Kiran Bharat Gaikwad
#Purpose: Estimate differentially expressed genes from rna-seq datasets of patients with FECD (DESeq2)

#Differential expression analysis
library("DESeq2")
file1 = read.csv("GSE142538_101724.txt", sep='\t', header = TRUE, row.names = 1, check.names =TRUE)
dim(file1)
head(file1)
countdata1 = as.matrix(file1)
dim(countdata1)
(colnames(countdata1))
condition1 = factor(c(rep("control",8),rep("fecd",14)))
condition1
coldata1 = data.frame(row.names = colnames(countdata1), condition1)
coldata1
ddsFull = DESeqDataSetFromMatrix(countData=countdata1, colData= coldata1, design =~ condition1)
ddsFull
dds = DESeq(ddsFull)
dds
res = results(dds)
res
summary(res)
res_ordered = res[order(res$pvalue),]
res_d = as.data.frame(res)
res_d
res_d$diffexpressed <- "NO"


#Upregulated and downregulated identification
#if log2FoldChange > 1 and padj < 0.05, set as "UP"
res_d$diffexpressed[res_d$log2FoldChange > 1.0 & res_d$padj < 0.05] <- "UP"
#if log2FoldChange < -1 and padj < 0.05, set as "DOWN"
res_d$diffexpressed[res_d$log2FoldChange < -1.0 & res_d$padj < 0.05] <- "DOWN"
res_d <- cbind(rownames(res_d), data.frame(res_d, row.names=NULL))
res_d
colnames(res_d)[1] <- "hgnc_symbol"
res_d

a = res_d[which(res_d$diffexpressed == 'UP'),]
b = res_d[which(res_d$diffexpressed == 'DOWN'),]


#Export results to files 
write.table(a, file='GSE142538_up.txt', sep='\t', row.names = FALSE, col.names=TRUE)
write.table(b, file='GSE142538_down.txt', sep='\t', row.names = FALSE, col.names=TRUE)
write.table(res_d, file='all_genes.tsv', sep='\t', row.names = FALSE, col.names = TRUE)