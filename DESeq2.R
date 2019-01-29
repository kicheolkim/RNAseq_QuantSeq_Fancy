library(readr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

### condition
which(condi$Sample == "JN08")
which(condi$Sample == "JN14")
which(condi$Sample == "JN24")
condi2 <- condi[-c(7,13,23),]

condi_P4_WT_CC <- condi2[condi2$condition=="P4_WT_CC",]
condi_P4_APC_CC <- condi2[condi2$condition=="P4_APC_CC",]
condi_P4_WT_SC <- condi2[condi2$condition=="P4_WT_SC",]
condi_P4_APC_SC <- condi2[condi2$condition=="P4_APC_SC",]

condi_P9_WT_CC <- condi2[condi2$condition=="P9_WT_CC",]
condi_P9_APC_CC <- condi2[condi2$condition=="P9_APC_CC",]
condi_P9_WT_SC <- condi2[condi2$condition=="P9_WT_SC",]
condi_P9_APC_SC <- condi2[condi2$condition=="P9_APC_SC",]


## prepare dataset
tmp_condi <- rbind(condi_P4_WT_SC, condi_P4_APC_SC, condi_P9_WT_SC, condi_P9_APC_SC)  # SC
tmp_counts <- jn_counts_all[,colnames(jn_counts_all) %in% tmp_condi$Sample]
row.names(tmp_counts) <- row.names(jn_counts_all)
tmp_condi$Sample; colnames(tmp_counts)

## DEG analysis
cds <- DESeqDataSetFromMatrix(countData = tmp_counts, colData = tmp_condi, design=~condition-1)
cds <- cds[ rowSums(counts(cds)) > 1, ]
cds <- DESeq(cds)


##### PCA plot
vst.cds <- vst(cds, blind=TRUE)

pcaMatrix <- assay(vst(cds, blind=TRUE))
pcaNorm <- prcomp(t(pcaMatrix), scale. = TRUE)
barplot(((pcaNorm$sdev)^2 / sum((pcaNorm$sdev)^2) * 100)[1:10], main = "Variance plot, scaled PCA for all samples")

plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(6, "Dark2")[condi_SC$condition], cex.axis=0.8, cex=1.3,
     pch = 16, xlab="PC1 (34.1%)", ylab="PC2 (16.4%)")
legend("bottomleft", c("P4 APC","P4 WT","P9 APC","P9 WT"), text.col = RColorBrewer::brewer.pal(6, "Dark2"), 
       col=RColorBrewer::brewer.pal(6, "Dark2"), pch=16, inset=.01, bty='n', cex=0.9)


##### Results for differential expressed gene analysis 
## results -1
results <- results(cds, alpha=0.05)
sum(results$padj < 0.05, na.rm=TRUE)
summary(results)

## results -2
res1 <- results(cds, alpha=0.05, contrast=c("condition","P4_APC_SC","P4_WT_SC"))   # P4 APC-WT
sum(res1$padj < 0.05, na.rm=TRUE)
summary(res1)
res1 <- as.data.frame(res1)
res1$GENEID <- row.names(res1)
res1 <- merge(res1, gene_list, by.x="GENEID", by.y="GENEID", sort=FALSE)

res2 <- results(cds, alpha=0.05, contrast=c("condition","P9_APC_SC","P9_WT_SC"))   # P9 APC-WT
sum(res2$padj < 0.05, na.rm=TRUE)
summary(res2)
res2 <- as.data.frame(res2)
res2$GENEID <- row.names(res2)
res2 <- merge(res2, gene_list, by.x="GENEID", by.y="GENEID", sort=FALSE)

res_both <- cbind(res1, res2)



### heatmap for significant genes
tmp_counts_plot <- assay(vst(cds))

# for heatmap
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = -sv, agglo.FUN=mean)
  as.hclust(dend)
}
# color for colorbar of heatmap
ann_colors = list(
  Type = c(WT=brewer.pal(12,"Paired")[9], APC=brewer.pal(12,"Paired")[10]),
  Time = c(P4= brewer.pal(12,"Paired")[11], P9=brewer.pal(12,"Paired")[12]))


# P4
res <- res1
tmp_plot_list <- subset(res, res$padj < 0.05); nrow(tmp_plot_list)  # & (res$log2FoldChange > 1 | res$log2FoldChange < -1)
tmp_plot_list <- tmp_plot_list[order(tmp_plot_list$padj),]
tmp_plot_list <- tmp_plot_list[1:200,]; nrow(tmp_plot_list)
tmp_sample_list <- subset(tmp_condi, time=="P4")
tmp_plot <- tmp_counts_plot[row.names(tmp_counts_plot) %in% tmp_plot_list$GENEID, colnames(tmp_counts_plot) %in% tmp_sample_list$Sample]
tmp_plot <- tmp_plot - rowMeans(tmp_plot)
tmp_plot_group <- data.frame(Type=tmp_sample_list$type, Time=tmp_sample_list$time)
tmp <- data.frame(geneid=row.names(tmp_plot))
tmp <- merge(tmp, tmp_plot_list, by.x="geneid", by.y="GENEID", sort=FALSE)
row.names(tmp_plot) <- make.names(tmp$GeneSymbol, unique=TRUE)
row.names(tmp_plot_group) <- colnames(tmp_plot)

pheatmap(tmp_plot, color = colorRampPalette(rev(brewer.pal(10,"RdBu")))(256), scale="row",
         fontsize=8, fontsize_row = 3.4, fontsize_col = 8, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         show_rownames=TRUE, show_colnames=FALSE, cluster_rows=TRUE, cluster_cols=FALSE, 
         cutree_rows=3, cutree_cols=1, clustering_callback=callback,
         annotation_col=tmp_plot_group, annotation_colors=ann_colors)


# P9
res <- res2
tmp_plot_list <- subset(res, res$padj < 0.05); nrow(tmp_plot_list) #& (res$log2FoldChange > 1 | res$log2FoldChange < -1)
tmp_plot_list <- tmp_plot_list[order(tmp_plot_list$padj),]
tmp_plot_list <- tmp_plot_list[1:200,]; nrow(tmp_plot_list)
tmp_sample_list <- subset(tmp_condi, time=="P9")
tmp_plot <- tmp_counts_plot[row.names(tmp_counts_plot) %in% tmp_plot_list$GENEID, colnames(tmp_counts_plot) %in% tmp_sample_list$Sample]
tmp_plot <- tmp_plot - rowMeans(tmp_plot)
tmp_plot_group <- data.frame(Type=tmp_sample_list$type, Time=tmp_sample_list$time)
tmp <- data.frame(geneid=row.names(tmp_plot))
tmp <- merge(tmp, tmp_plot_list, by.x="geneid", by.y="GENEID", sort=FALSE)
row.names(tmp_plot) <- make.names(tmp$GeneSymbol, unique=TRUE)
row.names(tmp_plot_group) <- colnames(tmp_plot)

pheatmap(tmp_plot, color = colorRampPalette(rev(brewer.pal(10,"RdBu")))(256), scale="none",
         fontsize=8, fontsize_row = 3.4, fontsize_col = 8, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         show_rownames=TRUE, show_colnames=FALSE, cluster_rows=TRUE, cluster_cols=FALSE, 
         cutree_rows=2, cutree_cols=1, clustering_callback=callback,
         annotation_col=tmp_plot_group, annotation_colors=ann_colors)



ann_colors = list(
  Type = c(WT=brewer.pal(12,"Paired")[9], APC=brewer.pal(12,"Paired")[10]),
  Time = c(P4= brewer.pal(12,"Paired")[11], P9=brewer.pal(12,"Paired")[12]))

# P4 & P9
tmp_plot_list <- subset(res, res$padj < 0.01)
tmp_plot <- tmp_counts_plot[row.names(tmp_counts_plot) %in% tmp_plot_list$GENEID,]
tmp_plot <- tmp_plot - rowMeans(tmp_plot)
tmp_plot <- tmp_plot[c("ENSMUSG00000020218.11","ENSMUSG00000053279.6","ENSMUSG00000036480.9","ENSMUSG00000020061.17","ENSMUSG00000020090.7","ENSMUSG00000033227.7","ENSMUSG00000000094.12","ENSMUSG00000007872.3"),]
tmp_plot_group <- data.frame(Type=tmp_condi$type, Time=tmp_condi$time)
tmp <- data.frame(geneid=row.names(tmp_plot))
tmp <- merge(tmp, tmp_plot_list, by.x="geneid", by.y="GENEID", sort=FALSE)
row.names(tmp_plot) <- make.names(tmp$GeneSymbol, unique=TRUE)
row.names(tmp_plot_group) <- colnames(tmp_plot)


pheatmap(tmp_plot, color = colorRampPalette(rev(brewer.pal(10,"RdBu")))(256), scale="none",
         fontsize=8, fontsize_row = 12, fontsize_col = 8, #breaks=seq(-3,4, by=0.0275), 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         show_rownames=TRUE, show_colnames=FALSE, cluster_rows=TRUE, cluster_cols=FALSE, 
         cutree_rows=2, cutree_cols=1, clustering_callback=callback,
         annotation_col=tmp_plot_group, annotation_colors=ann_colors)



#####################################################################################################################
# Spinal cord:     condi_P4_APC_SC, condi_P4_WT_SC, condi_P9_APC_SC, condi_P9_WT_SC

## prepare dataset
tmp_condi <- rbind(condi_P9_WT_SC, condi_P9_APC_SC)                              #################
tmp_counts <- jn_counts_all[,colnames(jn_counts_all) %in% tmp_condi$Sample]
row.names(tmp_counts) <- row.names(jn_counts_all)
tmp_condi$Sample; colnames(tmp_counts)

## DEG analysis
cds <- DESeqDataSetFromMatrix(countData = tmp_counts, colData = tmp_condi, design=~type)
cds <- cds[ rowSums(counts(cds)) > 1, ]
cds <- DESeq(cds)

results <- results(cds, alpha=0.05)
sum(results$padj < 0.05, na.rm=TRUE)
summary(results)

results <- as.data.frame(results)
results$GENEID <- row.names(results)
results <- merge(results, gene_list, by.x="GENEID", by.y="GENEID", sort=FALSE)

results <- results[order(results$padj),]
write.csv(results, "res_P9WT_CC-SC.csv")                                     #################


## heatmap for significant genes
tmp_counts_plot <- assay(vst(cds))
tmp_plot_list <- subset(results, results$padj < 0.05)
tmp_plot <- tmp_counts_plot[row.names(tmp_counts_plot) %in% tmp_plot_list$GENEID,]
tmp_plot <- tmp_plot - rowMeans(tmp_plot)
tmp_plot_group <- data.frame(time=tmp_condi$time, type=tmp_condi$type, tissue=tmp_condi$tissue)
tmp <- data.frame(geneid=row.names(tmp_plot))
tmp <- merge(tmp, tmp_plot_list, by.x="geneid", by.y="GENEID", sort=FALSE)
row.names(tmp_plot) <- make.names(tmp$GeneSymbol, unique=TRUE)
row.names(tmp_plot_group) <- colnames(tmp_plot)
pheatmap(tmp_plot, color = colorRampPalette(rev(brewer.pal(8,"RdBu")))(256), scale="none",
         fontsize=10, fontsize_row = 0.5, fontsize_col = 8, main="Heatmap for significant genes (adj.p < 0.05)", 
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, 
         annotation_col=tmp_plot_group)



