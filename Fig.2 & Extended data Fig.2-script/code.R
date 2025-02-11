#script for figure2
#..........................................plot correlation between spatial distance and genetic similarity..........................................
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(ggbeeswarm)
library(ggsignif)
library(patchwork)
library(ggpmisc)
setwd('./Desktop/R-data')
theme_set(ggpubr::theme_pubr()+theme(legend.position = "top"))
data<-read.table('data.txt',header=T,row.names = 1,check.names=FALSE);head(data)

chromosome<-colnames(data)[4:ncol(data)]
color_list<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C", "#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6")

pdf('./spearman_cor.pdf')
i=3
for (chro in chromosome){
  i=i+1
  p<-ggplot(data, aes(x = data$spatial_distance, y = data[,i]))+ geom_point(aes(color = patient)) +geom_rug(aes(color =patient)) +geom_smooth(aes(color = patient), method = lm,  se = FALSE, fullrange = TRUE)+scale_color_manual(values = color_list)+ggpubr::stat_cor(method='spearman',aes(color = patient))+labs(x='genetic_distance',y=chro)
  print(p)
}
dev.off()



#..........................................partial correlation..........................................
#install.packages("ppcor")
library(ppcor)
patient<-c('P01','P02','P03','P04','P05','P06','P09','P07')
data<-read.table('input-data.txt',header=T,row.names = 1);head(data)

pathway<-colnames(data)[4:ncol(data)]
i=3
j=0
string=character(ncol(data)*10)
for (pathway_id in pathway){
  print(pathway_id)
  i=i+1
  j=j+1
  pcor_total<-pcor.test(data$SNX20,
                        data[,i],
                        data[,c('tumor_purity')],
                        method="spearman")
  string[j]=paste(pathway_id,'all patients',pcor_total$p.value,pcor_total$estimate,sep=",")
  #print(paste('patient','pcor.pvalue','pcor.corr',sep=' '))
  for (key in patient){
    #print(key)
    j=j+1
    test1<-subset(data, `patient` ==key)
    pcor<-pcor.test(test1$SNX20,
                    test1[,i],
                    test1[,c('tumor_purity')],
                    method="spearman")
    string[j]=paste(pathway_id,key,pcor$p.value,pcor$estimate,sep=",")
  }
}
write.csv(file='resul.csv',paste(string))


#..........................................NMF..........................................
#step1
library(NMF)
library(BimodalIndex)
library(magrittr)
options(stringsAsFactors=F)

# arguments
args <- commandArgs(T)
cutoff <- args[1] %>% as.numeric()
med <- args[2]
top_n <- args[3] %>% as.numeric()
pfx <- ifelse(med == 'BiMeanSd', paste0('TNBC_', cutoff, 'pct-tumor_', med),
        paste0('TNBC_', cutoff, 'pct-tumor_', med, '-top', top_n)
)
output <- paste0('output_TNBC_', cutoff, 'pct-tumor_', med)
dir.create(output)
# prepare input data
## expression data
tnbc_all <- readRDS('./TNBC_expr_TPM.Rds')
str(tnbc_all)
rownames(tnbc_all) <- tnbc_all$GeneSymbol
tnbc_all$GeneSymbol <- NULL
str(tnbc_all)
## expr_gene_num > 5000 TNBC cel clusters
gene_num <- apply(tnbc_all, 2, function(x) { length(x[x>0]) })
str(gene_num)
gt_5k_sam <- gene_num[gene_num > 5000] %>% names()
str(gt_5k_sam)
## subset samples based on tumor purity
meta_data <- read.delim('./TNBC_metadata.tsv',
        header = T, check.names = F, row.names = 4) %>% .[gt_5k_sam, ]
str(meta_data)
tumor_samples <- subset(meta_data, `Tumor (%)` > cutoff) %>% rownames
tnbc_tumor <- tnbc_all[, tumor_samples] %>% as.matrix
## gene_qc
num_cell <- apply(tnbc_tumor, 1, function(x) { length(x[x > 0]) })
ge_3_genes <- num_cell[num_cell >= 3] %>% names()
str(ge_3_genes)
tnbc_tumor <- tnbc_tumor[ge_3_genes, ]
str(tnbc_tumor)
## extract top_n genes
BiMeanSd_fun <- function(x) {
  bi <- bimodalIndex_cj(x)
  bi$mean <- apply(x, 1, mean)
  bi$sd <- apply(x, 1, sd)
  return(bi)
}
### calculate & select top_n genes
gene_val <- switch(med,
  'cv' = apply(tnbc_tumor, 1, sd) / rowMeans(tnbc_tumor),
  'mad' = apply(tnbc_tumor, 1, mad),
  'sd' = apply(tnbc_tumor, 1, sd),
  'BiMeanSd' = BiMeanSd_fun(tnbc_tumor)
)
gene_sel <- 'init'
if (med == 'BiMeanSd') {
  gene_sel <- subset(gene_val, BI >= 1.5 & mean >= quantile(gene_val$mean, .25) & sd >= quantile(gene_val$sd, .5)) %>%
        rownames()
} else {
  gene_sel <- sort(gene_val, decreasing=T)[1:top_n] %>% names()
}
expr <- tnbc_tumor[gene_sel,, drop=F ]
str(expr)
# Estimation of the rank
estim.r <- nmf(expr, 2:8, nrun = 50, seed = 123456, .opt = "vp25")
pdf(paste0(output, '/', pfx, '_nmf_estimate-rank.pdf'), width = 15, height = 15)
plot(estim.r)
consensusmap(estim.r)
dev.off()
save(expr, gene_sel, gene_val, estim.r, meta_data, file = paste0(output, '/', pfx, '_nmf_estimate-rank.Robj'))

#step2
# arguments
args <- commandArgs(T)
cutoff <- args[1] %>% as.numeric()
med <- args[2]
top_n <- args[3] %>% as.numeric()
k <- args[4] %>% as.numeric()
pfx <- ifelse(med == 'BiMeanSd', paste0('TNBC_', cutoff, 'pct-tumor_', med),
     paste0('TNBC_', cutoff, 'pct-tumor_', med, '-top', top_n)
)
output <- paste0('output_TNBC_', cutoff, 'pct-tumor_', med)
load(paste0(output, '/', pfx, '_nmf_estimate-rank.Robj'))
# data
## annotation
ann_col <- meta_data[colnames(expr), c('Patient', 'Site', 'Tumor (%)', 'Immue cell (%)', 'Fibroblast (%)')]
# Run NMF based on selected K
nmf_res <- nmf(expr, rank = k, nrun=200, seed=123456, .opt="vp25")
# plot
pdf(paste0(output, '/', pfx, '_nmf_rank', k, '_basis-coef-consensus_plot.pdf'), width=18, height=8)
layout(matrix(c(1,2,3,3), nrow=1, byrow=F))
# basis components  labRow=NA
basismap(nmf_res)
# mixture coefficients
coefmap(nmf_res, labCol=NA)
# consensus plotting
conses_map <- consensusmap(nmf_res, annCol = ann_col, labRow=NA, labCol=NA)
dev.off()
save(expr, nmf_res, meta_data, conses_map, ann_col, file = paste0(output, '/', pfx, '_nmf_rank', k, '_result.Robj'))

#step3
library(NMF)
library(RColorBrewer)
library(dendextend)
library(magrittr)
library(pheatmap)
options(stringsAsFactors=F)
rm(list=ls())
k <- 5
cutoff <- 75
med <- 'sd'
top_n <- 2000
n_gene <- 100
setwd('./output_TNBC_75pct-tumor_sd')
pfx <- paste0('TNBC_', cutoff, 'pct-tumor_', med, '-top', top_n)

# using built-function: extractFeatures
load(paste0(pfx, '_nmf_rank', k, '_result.Robj'))
s <- extractFeatures(nmf_res)
lapply(s, function(x) { rownames(expr)[x] })

# scale method: make rowSum to 1
gene_sel <- basis(nmf_res) %>%
        NMF:::scale_mat(., scale = 'r1') %>%
        apply(., 2, function(x) {
        sort(x, decreasing=T) %>% head(n_gene) %>% names()
}) %>% as.character %>% unique
# get annotation
## get cluster ID
cluster <- predict(nmf_res)
cluster_info <- data.frame(row.names = names(cluster), Cluster = as.character(cluster))
ann_col <- cbind(ann_col[c(1:3)], cluster_info[rownames(ann_col), 'Cluster', drop=F])
## annotation colors
ann_color <- list(
  Patient = setNames(brewer.pal(12, "Paired")[1:8], sort(unique(ann_col$Patient))),
  Site = setNames(c('#3d3d3d', '#919191', 'lightgrey'), c('metastasis-D', 'metastasis-A', 'primary')),
  Cluster = setNames(brewer.pal(9, "Set1")[1:k], c(1:k)),
  `Tumor (%)` = setNames(colorRampPalette(c("#DEEBF7", "blue"))(101), 0:100)

)
## load expression data
expr_scale <- log2(expr + 1) %>% t() %>% scale %>% t()
## colors
cutree_k <- 5
bk <- seq(-2, 2, 0.5)
my_cols <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))(length(bk))
col_order <- rownames(ann_col[order(ann_col$Cluster), ])
pdf(paste0( pfx, '_nmf_rank', k, '_nmf-selected-top', n_gene, '-genes_scale-r1_samOrderByCls.pdf'), 12, 8.5)
pm <- pheatmap(expr_scale[gene_sel, col_order], border_color = NA,
    breaks = bk, color = my_cols,
    annotation_col = ann_col, annotation_colors = ann_color,
    show_rownames = F, show_colnames = F,
    cluster_rows = T, cluster_cols = F,
    clustering_method = "ward.D2",
    cutree_rows = cutree_k
)
dev.off()
# color branches and extract leaves
dend <- pm$tree_row %>%
    as.dendrogram() %>%
    set("labels_to_character") %>%
    color_branches(k = cutree_k, groupLabels = TRUE)
# Prune cutree to dendlist
clusters <- dendextend::cutree(dend, cutree_k, order_clusters_as_data = FALSE)
dendl <- vector("list", cutree_k)
for (i in 1:cutree_k) {
  leves_to_prune <- labels(dend)[clusters != i]
  dendl[[i]] <- prune(dend, leves_to_prune)
}
names(dendl) <- paste0('Module', 1:length(dendl))
gene_module <- sapply(dendl, labels)
str(gene_module)

save(pm, dendl, gene_module, file = paste0( pfx, '_nmf_rank', k, '_nmf-selected-top', n_gene, '-genes_scale-r1_module.Robj'))
all_line <- c()
for (n in names(gene_module)) {
  line <- paste(c(n, gene_module[[n]]), collapse = '\t')
  all_line <- c(all_line, line)
}
writeLines(all_line, con = paste0( pfx, '_nmf_rank', k, '_nmf-selected-top', n_gene, '-genes_scale-r1_module.txt'))
## plot dendrogram
pdf(paste0( pfx, '_nmf_rank', k, '_nmf-selected-top', n_gene, '-genes_scale-r1_dendrogram.pdf'), 4, 10)
plot(rev(dend), main = "Row dendrogram", horiz=T)
dev.off()
# enrichment ----------------
# module enrichment analysis ----------------
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
options(stringsAsFactors=F)
## database
H_signature <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(H_signature)

res_hallmark <- lapply(names(gene_module), function(x) {
  res <- enricher(gene_module[[x]], TERM2GENE=H_signature) %>% as.data.frame
  if ( nrow(res) > 0 ) res$Group <- x
  return(res)
}) %>% do.call(rbind, .)
head(res_hallmark)
### plot
pdf(paste0(pfx, '_nmf_rank', k, '_nmf-selected-top', n_gene, '-genes_scale-r1_module_enrichHALLMARK.pdf'), 6, 4)
ggplot(res_hallmark, aes(Group, ID)) + geom_point(aes(color = p.adjust, size = Count)) +
    scale_colour_distiller(palette = 'Reds') +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color='black'),
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.ticks = element_line(color = 'black'),
          axis.title = element_blank())

dev.off()
## GO enrichment ------------------
res_go  <- lapply(names(gene_module), function(x) {
    print(paste('Running for', x))
    res <- enrichGO(gene  = gene_module[[x]],
          keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
          ont = "BP", pAdjustMethod = "BH",
          pvalueCutoff  = 0.01, qvalueCutoff  = 0.05
    ) %>% as.data.frame()

    if ( nrow(res) > 0 ) res$Group <- x
    return(res)
})  %>% do.call(rbind, .)
pdf(paste0( pfx, '_nmf_rank', k, '_nmf-selected-top', n_gene, '-genes_scale-r1_module_enrichGO.pdf'), 8, 19)
ggplot(res_go, aes(Group, Description)) + geom_point(aes(color = p.adjust, size = Count)) +
    scale_colour_distiller(palette = 'Reds') +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color='black'),
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.ticks = element_line(color = 'black'),
          axis.title = element_blank())

dev.off()
## KEGG metabolism enrichment -------------
### metabolism pathway
f_metabolism <- readLines('../input/KEGG_metabolism_pathway.tsv')
pw_df <- lapply(f_metabolism, function(l) {
    ele <- strsplit(l, '\t')[[1]]
    data.frame(gs_name = rep(ele[1], length(ele) - 2),
        gene_symbol = ele[-c(1:2)]
    )
}) %>% do.call(rbind, .)
res_metabolism <- lapply(names(gene_module), function(x) {
  res <- enricher(gene_module[[x]], TERM2GENE=pw_df) %>% as.data.frame
  if ( nrow(res) > 0 ) res$Group <- x
  return(res)
}) %>% do.call(rbind, .)
head(res_metabolism)
### plot
pdf(paste0(pfx, '_nmf_rank', k, '_nmf-selected-top', n_gene, '-genes_scale-r1_module_enrichMetabolism.pdf'), 5, 3.5)
ggplot(res_metabolism, aes(Group, ID)) + geom_point(aes(color = p.adjust, size = Count)) +
    scale_colour_distiller(palette = 'Reds') +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color='black'),
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.ticks = element_line(color = 'black'),
          axis.title = element_blank())

dev.off()
save(res_hallmark, res_go, res_metabolism, file = paste0(pfx, '_nmf_rank', k, '_nmf-selected-top', n_gene, '-genes_scale-r1_module_enrichment.Robj'))
