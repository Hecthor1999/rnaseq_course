#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)

#load feature counts table and metadata
counts_data <- read.table("featurecounts.txt", header= TRUE, row.names = 1)
col_data <- read.csv("metadata.txt", header = TRUE, row.names = 1)

# List of columns to exclude
columns_to_exclude <- c("Chr", "Start", "End", "Strand", "Length")
# Subset to keep only the relevant columns for DESeq2
counts_data_subset <- counts_data[, !(colnames(counts_data) %in% columns_to_exclude)]
# Create a vector of new column names 
new_column_names  <- c("HER21", "HER22", "HER23", "NonTNBC1", "NonTNBC2", "NonTNBC3","Normal1","Normal2","Normal3","TNBC1","TNBC2", "TNBC3")
# Assign the new row names to counts_data
colnames(counts_data_subset) <- new_column_names

#Make sure col_data has the same names as counts_data_subset
all(colnames(counts_data_subset) %in% rownames(col_data))
# and make sure they are in the same order
all(colnames(counts_data_subset) == rownames(col_data))

#construct a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_data_subset, colData = col_data, design = ~ subtypes)
dds

#set the reference level
dds$subtypes <- relevel(dds$subtypes, ref = "Normal")

#run Deseq2
dds <- DESeq(dds)
results <- results(dds)
results

#explore results
summary(results)
resultsNames(dds)

#Plot counts
plotCounts(dds, gene=which.min(results$padj), intgroup="subtypes", pch=19)

#extract transformed values to remove dependence of the variance on the mean
vsd <- vst(dds, blind=TRUE)
head(assay(vsd), 3)

# visualization
ntd <- normTransform(dds)
#BiocManager::install("vsn")
library(vsn)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

#PCA plots
plotPCA(vsd, intgroup=c("cancer", "subtypes"))


#contrast
res_1 <- results(dds, contrast=c("subtypes", "TNBC", "HER2"), alpha = 0.05) 
summary(res_1) 
sum( res_1$padj < 0.05, na.rm=TRUE) # number of differentially expressed genes
res_2 <- results(dds, contrast=c("subtypes", "TNBC", "NonTNBC"), alpha = 0.05)
summary(res_2)
sum( res_2$padj < 0.05, na.rm=TRUE)
res_3 <- results(dds, contrast=c("subtypes", "NonTNBC", "HER2"), alpha = 0.05)
summary(res_3)
sum( res_3$padj < 0.05, na.rm=TRUE)

#get the gene ids
#1
library( "biomaRt" )
res_1$ensembl <- rownames(res_1)
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res_1$ensembl,
                  mart = ensembl )
idx <- match( res_1$ensembl, genemap$ensembl_gene_id)
res_1$entrez <- genemap$entrezgene_id[ idx ]
res_1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

#2
res_2$ensembl <- rownames(res_2)
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res_2$ensembl,
                  mart = ensembl )
idx <- match( res_2$ensembl, genemap$ensembl_gene_id)
res_2$entrez <- genemap$entrezgene_id[ idx ]
res_2$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

#3
res_3$ensembl <- rownames(res_3)
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res_3$ensembl,
                  mart = ensembl )
idx <- match( res_3$ensembl, genemap$ensembl_gene_id)
res_3$entrez <- genemap$entrezgene_id[ idx ]
res_3$hgnc_symbol <- genemap$hgnc_symbol[ idx ]


#create a subset of the results with the differentially expressed genes
res_TNBC_HER2 <- res_1[ which(res_1$padj < 0.05), ]
res_TNBC_NonTNBC <- res_2[ which(res_2$padj < 0.05), ]
res_NonTNBC_HER2 <- res_3[ which(res_3$padj < 0.05), ]

#enhanced volcano plots

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

EnhancedVolcano(res_TNBC_HER2,
                lab = res_TNBC_HER2$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'TNBC vs HER2',
                pointSize = 2.0,
                labSize = 3.0,
                colAlpha = 2,
                legendPosition = 'None',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

EnhancedVolcano(res_TNBC_NonTNBC,
                lab = res_TNBC_NonTNBC$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'TNBC vs NonTNBC',
                pointSize = 2.0,
                labSize = 3.0,
                colAlpha = 2,
                legendPosition = 'None',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

EnhancedVolcano(res_NonTNBC_HER2,
                lab = res_NonTNBC_HER2$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'NonTNBC vs HER2',
                pointSize = 2.0,
                labSize = 3.0,
                colAlpha = 2,
                legendPosition = 'None',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

#expression of individual genes

CDKN2A <- plotCounts(dds, "ENSG00000147889", intgroup = c("subtypes"), returnData = TRUE)
boxplot(count ~ subtypes, data=CDKN2A, main="Expression of CDKN2A", col="skyblue")

#calculate the mean and standard deviation of the read counts in TNBC
mean_CDKN2A_TNBC <- mean(CDKN2A[10:12,1])
sd_CDKN2A_TNBC <- sd(CDKN2A[10:12,1])

#calculate the mean and standard deviation of the of the read counts in rest of the samples combined
mean_CDKN2A_HER2_NonTNBC_Normal <- mean(CDKN2A[1:9,1])
sd_CDKN2A_HER2_NonTNBC_Normal <- sd(CDKN2A[1:9,1])

SNORA53 <- plotCounts(dds, "ENSG00000212443", intgroup = c("subtypes"), returnData = TRUE)
boxplot(count ~ subtypes, data=SNORA53, main="Expression of SNORA53", col="skyblue")

#calculate the mean and standard deviation of the read counts in HER2
mean_SNORA53_TNBC_NonTNBC_Normal <- mean(SNORA53[1:3,1])
sd_SNORA53_TNBC_NonTNBC_Normal <- sd(SNORA53[1:3,1])

#calculate the mean and standard deviation of the of the read counts in rest of the samples combined
mean_SNORA53_HER2 <- mean(SNORA53[4:12,1])
sd_SNORA53_HER2 <- sd(SNORA53[4:12,1])

ESR1 <- plotCounts(dds, "ENSG00000091831", intgroup = c("subtypes"), returnData = TRUE)
boxplot(count ~ subtypes, data=ESR1, main="Expression of ESR1", col="skyblue")

#calculate the mean and standard deviation of the read counts in NonTNBC
mean_ESR1_NonTNBC <- mean(ESR1[4:6,1])
sd_ESR1_NonTNBC <- sd(ESR1[4:6,1])

#calculate the mean and standard deviation of the of the read counts in rest of the samples combined
mean_ESR1_TNBC_HER2_Normal <- mean(ESR1[c(1:3, 7:12),1])
sd_ESR1_TNBC_HER2_Normal <- sd(ESR1[c(1:3, 7:12),1])

#export the results:
write.csv( as.data.frame(res_1), file="results_TNBC_vs_HER2.csv")
write.csv( as.data.frame(res_2), file="results_TNBC_vs_NonTNBC.csv")
write.csv( as.data.frame(res_3), file="results_NonTNBC_vs_HER2.csv")

write.csv( as.data.frame(res_TNBC_HER2), file="results_TNBC_vs_HER2_filtered.csv")
write.csv( as.data.frame(res_TNBC_NonTNBC), file="results_TNBC_vs_NonTNBC_filtered.csv")
write.csv( as.data.frame(res_NonTNBC_HER2), file="results_NonTNBC_vs_HER2_filtered.csv")


#BiocManager::install("AnnotationDbi")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


universe1 <- rownames(res_1) #all of the genes
genes1 <- rownames(res_TNBC_HER2) #differentially expressed genes

ego1 <- enrichGO(gene          = genes1,
                 universe      = universe1,
                 OrgDb         = org.Hs.eg.db,
                 keyType = "ENSEMBL",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)


universe2 <- rownames(res_2)
genes2 <- rownames(res_TNBC_NonTNBC)

ego2 <- enrichGO(gene          = genes2,
                 universe      = universe2,
                 OrgDb         = org.Hs.eg.db,
                 keyType = "ENSEMBL",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

universe3 <- rownames(res_3)
genes3 <- rownames(res_NonTNBC_HER2)

ego3 <- enrichGO(gene          = genes3,
                 universe      = universe3,
                 OrgDb         = org.Hs.eg.db,
                 keyType = "ENSEMBL",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)


write.csv(as.data.frame(ego1), file="Results_enrichGO_TNBC_vs_HER2.csv")
write.csv(as.data.frame(ego2), file="Results_enrichGO_TNBC_vs_NonTNBC.csv")
write.csv(as.data.frame(ego3), file="Results_enrichGO_NonTNBC_vs_HER2.csv")

#BiocManager::install("enrichplot")
library(enrichplot)

# Only included plots for contrast TNBC vs HER2 due to space constraints 
dotplot(ego1, showCategory=10) + ggtitle("TNBC vs HER2")
upsetplot(ego1) + ggtitle("TNBC vs HER2")

