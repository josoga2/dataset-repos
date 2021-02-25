#---------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(dendextend)
#install.packages("colorRamp2")
library(colorRamps)

library(DESeq2)


#library(GeoTcgaData)

#BiocManager::install("GeoTcgaData")


#rawdData <- read.csv('https://github.com/wale-joseph/single-cell-notebooks/raw/main/drosophila_counted.csv', sep='\t')

rawdData <- read.csv('C:/Users/Wale Joseph/Downloads/deseqFile.csv', sep = ',')
head(rawdData)


#data may potentially contains non-numbers 
#rawdData <- rawdData[!is.na(as.numeric(as.character(rawdData$panc_rep1, rawdData$panc_rep2,
#                                                    rawdData$neur_rep1, rawdData$neur_rep2))),]

metadata <-read.table('C:/Users/Wale Joseph/Downloads/design.txt', sep='\t', header = T )
#Construct DESEQ2 object
dds <- DESeqDataSetFromMatrix(countData = rawdData, colData = metadata, design = ~dex, tidy=T)

head(dds)

deseqOutput <- DESeq(dds)

outResult <- results(deseqOutput)

head(outResult) 

#export datasets
upReg <- subset(outResult, log2FoldChange >2, padj <= 0.05) #You can increase the lfc, this is a comparison between 2 cell types.  
downReg <- subset(outResult, log2FoldChange <2, padj <= 0.05)

write.table(upReg,  file = 'upReg.tsv', sep = '\t')
write.table(downReg,  file = 'downReg.tsv', sep = '\t')


#visualizations 
#heatmap  ----you can reduce the amount of genes, if you are on a personal PC 
myMatrix <- as.matrix(rawdData[ ,c(2:5)])
myMatrix <- t(myMatrix)
class(myMatrix)

dend = hclust(dist(myMatrix,method="maximum"),method="ward.D")

par(mfrow=c(1,1))
Heatmap(t(myMatrix), cluster_columns = T,
        row_names_side = "left",
        row_dend_side = "left",
        row_names_gp=gpar(cex=1),
        row_dend_width = unit(1, "cm"),
        clustering_distance_rows ="kendall",
        clustering_method_rows = "complete",
        km=2, height = 1, heatmap_legend_param = list(title='LogFC'))




with(outResult, plot(log2FoldChange, -log10(pvalue), col='blue'))

par(mfrow=c(1,1))
with(subset(outResult, padj<0.05), points(log2FoldChange, -log10(padj), col='blue'))
with(subset(outResult, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj),pch=20, col='red', main="Volcano Plot"))
abline(v=2, col='red')
abline(v=-2, col='green')
with(subset(outResult, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj),pch=20, col='red', labels = c('gene1', 'gene2', 'gene3', 'gene4') ))
