#script for figure1
#..........................seurat clustering...........................
setwd('./')
source('tnbc_seurat_clustering-preset.r')
# ----------- step1 -------------
if(1 %in% steps){
print(paste('run step1 at', Sys.time(), '------'))
## load expr matrix
tpm <- readRDS('../TPM.Rds')
rownames(tpm) <- tpm$GeneSymbol
tpm$GeneSymbol <- NULL
dim(tpm)
### create object 
pbmc <- CreateSeuratObject(counts=as.matrix(tpm), project=pfx, min.cells=minCell, min.features=geneLow)
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, pattern='^MT-')
## set group info
pbmc$locus_ori <- sapply(rownames(pbmc@meta.data), function(x) { strsplit(x, '-')[[1]][2] })
pbmc$locus <- sub('A|D', 'metastasis', pbmc$locus_ori) %>% sub('BC', 'primary', .)
pbmc$patient <- sub('-.*', '', rownames(pbmc@meta.data))
## qc
pbmc <- subset(pbmc, subset=nFeature_RNA > geneLow & percent.mt < mitoHigh)
dim(pbmc)
str(pbmc@meta.data)
# visual features
pdf(paste0(outdir, pfx, '_step1_createObject.pdf'), 10)
print(VlnPlot(pbmc, features=c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol=3))
plot1 <- FeatureScatter(pbmc, feature1='nCount_RNA', feature2='percent.mt')
plot2 <- FeatureScatter(pbmc, feature1='nCount_RNA', feature2='nFeature_RNA')
CombinePlots(plots=list(plot1, plot2))
dev.off()
### data precessing
# normlize data
pbmc <- NormalizeData(pbmc, normalization.method='LogNormalize', scale.factor=10000)
save(pbmc, file=paste0(outdir, pfx, '_step1s_norm_data.Robj'))
print(paste('Done for step1 at', Sys.time(), '------'))
}
# ----------- step2 -----------
if(2 %in% steps){
print(paste('run step2 at', Sys.time(), '------'))
# Find variable gene
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot')
str(VariableFeatures(object=pbmc))
# Scaling the data
pbmc <- ScaleData(pbmc, features=rownames(pbmc))
# run PCA
pbmc <- RunPCA(pbmc, features=VariableFeatures(object = pbmc))
print(pbmc[['pca']], dims=1:5, nfeatures=5)
pdf(paste0(outdir, pfx, '_step2_runPCA.pdf'), 10)
print(DimPlot(pbmc, reduction = "pca", group.by = 'patient', cols = randomColor))
opar <- par(mfrow=c(2,3))
print(DimHeatmap(pbmc, dims=1:6, cells=200, balanced=TRUE))
print(DimHeatmap(pbmc, dims=7:12, cells=200, balanced=TRUE))
print(DimHeatmap(pbmc, dims=13:18, cells=200, balanced=TRUE))
print(DimHeatmap(pbmc, dims=19:24, cells=200, balanced=TRUE))
print(DimHeatmap(pbmc, dims=25:30, cells=200, balanced=TRUE))
print(ElbowPlot(pbmc, ndims=30))
dev.off()
save(pbmc, file=paste0(outdir, pfx, '_step2_runPCA.Robj'))
print(paste('Done for step2 at', Sys.time(), '------'))
}
# ----------- step3 --------------
if(3 %in% steps){
print(paste('run step3 at', Sys.time(), '------'))
sfx <- paste0('res', nReso, '_', nPcDim, 'pc')
load(paste0(outdir, pfx, '_step2_runPCA.Robj'))
### Cluster cells
pbmc <- FindNeighbors(pbmc, dims=1:as.numeric(nPcDim), k.param = 30)
pbmc <- FindClusters(pbmc, resolution=as.numeric(nReso))
# visual
pbmc <- RunTSNE(pbmc, dims=1:as.numeric(nPcDim))
pbmc <- RunUMAP(pbmc, dims=1:as.numeric(nPcDim))
nIdent <- nlevels(Idents(pbmc))
nSample <- length(unique(pbmc@meta.data$patient))
str(randomColor)
pdf(paste0(outdir, pfx, '_step3_cluster_', sfx, '.pdf'), 8.5)
print(DimPlot(pbmc, reduction='tsne', group.by='patient', cols=randomColor[1:nSample]))
print(DimPlot(pbmc, reduction='tsne', group.by='locus', cols=randomColor[1:2]))
print(DimPlot(pbmc, reduction='tsne', cols=randomColor, label=T))
dev.off()
data <- as.matrix(GetAssayData(pbmc))
meta <- pbmc@meta.data
tsne <- Embeddings(pbmc, reduction='tsne')
varGene <- VariableFeatures(object=pbmc)
save(data, meta, tsne, varGene, file=paste0(outdir, pfx, '_step3_data-meta-tsne-varGene_', sfx, '.Robj'))
save(pbmc, file=paste0(outdir, pfx, '_step3_cluster_', sfx, '.Robj'))
### plot known markers
allmarkers <- list(
  Epithelial = c('KRT19', 'EPCAM', 'KRT18', 'KRT14'),
  Immune = c('PTPRC', 'LAPTM5', 'IL2RG'),
  Fibroblast = c('FAP', 'COL1A1', 'COL3A1', 'DCN')
)
pdf(paste0(outdir, pfx, '_step3_plotKnownMarker_', sfx, '.pdf'), width = 10)
# feature and violin plot
for (markersCandidate in allmarkers)
{
  print('want to VlnPlot')
  print(markersCandidate)
  markers.toplot <- intersect(rownames(pbmc), markersCandidate)
  if (length(markers.toplot) > 0)
  {
    print(VlnPlot(pbmc, features = markers.toplot))
    print(FeaturePlot(pbmc, features = markers.toplot, reduction = 'tsne'))
  }
}
dev.off()
### find markers
pbmcMarkers <- FindAllMarkers(pbmc)
str(pbmcMarkers)
# save top10 top80 top25
top80 <- pbmcMarkers %>% group_by(cluster) %>% top_n(80, avg_logFC)
top10 <- pbmcMarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top25 <- pbmcMarkers %>% group_by(cluster) %>% top_n(25, avg_logFC)
top25 <- as.data.frame(top25)
top25 <- top25[order(top25$cluster, -top25$avg_logFC),]
top25Df <- data.frame(row.names=1:25)
for(c in unique(top25$cluster))
{
  genes <- top25[top25$cluster == c,]$gene
  top25Df[, paste0('c', c)] <- c(genes, rep('NA', 25-length(genes)))
}
write.csv(top25Df, paste0(outdir, pfx, '_step3_findMarkers_top25_', sfx, '.csv'), row.names=F)
write.csv(top80, paste0(outdir, pfx, '_step3_findMarkers_top80_', sfx, '.csv'), row.names=F)
save(pbmcMarkers, top10, top80, top25, file=paste0(outdir, pfx, '_step3_findMarkers_', sfx, '.Robj'))
# do heatmap
pbmc <- ScaleData(pbmc, features=c(VariableFeatures(pbmc), top10$gene))
pdf(paste0(outdir, pfx, '_step3_findMarkers_top10heatmap_', sfx, '.pdf'), width=1.5*nIdent, height=0.95*nIdent)
print(DoHeatmap(pbmc, features=top10$gene, lines.width=10) + NoLegend())
dev.off()

#system(paste0('sh /home/zhangcuijuan/tool/others/pdf2png.sh ', outdir, pfx, '_step3_plotKnownMarker_', sfx, '.pdf'))
print(paste('Done for step3 at', Sys.time(), '------'))
}

if(4 %in% steps){

sfx <- paste0('res', nReso, '_', nPcDim, 'pc')
load(paste0(outdir, pfx, '_step3_cluster_', sfx, '.Robj'))

cluster <- read.delim(paste0(outdir, 'cluster_res4.txt'), header=T, as.is=T)
str(cluster)

pbmc$celltype.labels <- plyr::mapvalues(x = Idents(pbmc), from=cluster$cluster, to=cluster$celltype)
#pbmc$newLabels <- 'Immune'
#pbmc@meta.data[pbmc@meta.data$celltype.labels == 'Epithelial cell', 'newLabels'] <- 'Epithelial'
#pbmc@meta.data[pbmc@meta.data$celltype.labels == 'CAF', 'newLabels'] <- 'Stromal'


nCelltype <- nlevels(pbmc@meta.data$celltype.labels)
nIdent <- nlevels(Idents(pbmc))
nSample <- nlevels(pbmc@meta.data$orig.ident)

pdf(paste0(outdir, pfx, '_step4_label_celltype_', sfx, '.pdf'), 8.5)
print(DimPlot(pbmc, reduction='tsne', group.by='orig.ident', cols=randomColor[1:nSample]))
print(DimPlot(pbmc, reduction='tsne', cols=randomColor[1:nIdent], label=T))
print(DimPlot(pbmc, reduction='tsne', group.by='celltype.labels', cols=randomColor[1:nCelltype]))

dev.off()

save(pbmc, file=paste0(outdir, pfx, '_step4_label_celltype_', sfx, '.Robj'))
}




#..........................deconstructSigs...........................
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
setwd('./')

library(BSgenome.Hsapiens.UCSC.hg38)
library(deconstructSigs)
data<-read.table('1-raw-data-hg38-Mut_equal_more_30.txt',header=T)
head(data)
class(data)
sigs.input <- mut.to.sigs.input(mut.ref = data, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg =BSgenome.Hsapiens.UCSC.hg38)
dim(sigs.input)
head(sigs.input)
write.csv(sigs.input,file='2-sig_input-Mut_equal_more_30.csv')

w=lapply(unique(rownames(sigs.input)) , function(i){
  sample_1 = whichSignatures(tumor.ref = sigs.input[,], 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'exome2genome')
  print(i)
  return(sample_1$weights)
})
w=do.call(rbind,w)
head(w)
write.csv(w,file='3-result-Mut_equal_more_30.csv')
library(pheatmap)
pdf('mutation-signature-cosmic-v2-Mut_equal_more_30.pdf')
pheatmap(t(w),cluster_rows = F,cluster_cols = T,show_colnames=F)
dev.off()