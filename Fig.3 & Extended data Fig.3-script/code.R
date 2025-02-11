#script for figure3

#.....................................copy number variation in 4 subtypes...........................
library(pheatmap)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggrepel)
library(dplyr)
library(tidyr)
library(reshape2)
#mut <- read.table("data_mutations.txt",header=T,skip = 1,sep="\t",stringsAsFactors=F)
mut <- read.table("../all_thresholded.by_genes.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F,row.names=1)
drive <- read.table("./Breast_cancer_driver_genes_f.txt",header=T,sep="\t",stringsAsFactors=F)
subtype <- read.table("../ann_col.csv",header=T,stringsAsFactors=F,row.names=1,sep=",")
rownames(subtype) <- gsub("-","_",rownames(subtype))
aa <- intersect(rownames(subtype),colnames(mut))
mut_f <- mut[intersect(drive$Gene,rownames(mut)),aa]
subtype <- subtype[aa,]
anno_col <- data.frame(Cluster=rownames(subtype))
#anno_col$Patient <- sub("_.*","",anno_col$Cluster)
rename <- read.table("./rename_match.txt",sep="\t",header=T)
anno_col$Subtype <- subtype$Cluster
anno_col$Subtype <- paste("Cluster",anno_col$Subtype,sep="")
anno_col$Subtype <- gsub("Cluster1","IM regulatory",anno_col$Subtype)
anno_col$Subtype <- gsub("Cluster2","OXPHOS",anno_col$Subtype)
anno_col$Subtype <- gsub("Cluster3","LAR",anno_col$Subtype)
anno_col$Subtype <- gsub("Cluster4","IM/MES",anno_col$Subtype)
anno_col$Patient <- sub("_.*","",anno_col$Cluster)
anno_col$Patient %<>% plyr::mapvalues(from = rename$Patient_name,to = rename$Patient)
anno_col$Tissue <- sapply(strsplit(as.character(anno_col$Cluste),"_",fixed=TRUE),"[",2)
anno_col$Tissue <- gsub("A","metastasis_A",anno_col$Tissue)
anno_col$Tissue <- gsub("D","metastasis_D",anno_col$Tissue)
anno_col$Tissue <- gsub("BC","primary",anno_col$Tissue)
rownames(anno_col) <- anno_col$Cluster
anno_col$Cluster <- NULL
anno_col <- anno_col %>% mutate(Subtype=fct_relevel(Subtype,c("IM regulatory","IM/MES","LAR","OXPHOS"))) %>% arrange(Subtype)
OX_sam <- anno_col[which(anno_col$Subtype == "OXPHOS"),]
mut_f_OX <- mut_f[,rownames(OX_sam)]
#f <- function(x) sum(x!=0)
f <- function(x) sum( x < 0 | x > 0)
count <- as.data.frame(apply(mut_f_OX,1,f))
colnames(count) <- "count"
count$gene <- rownames(count)
drive_f <- rownames(count[which(count[,1] > ceiling(ncol(mut_f_OX) * 0.5)),])
length(drive_f)
pat_mut_1 <- mut_f[drive_f,rownames(anno_col)]
#anno_col$Subtype <- anno_col$group
#gene amp p-value
allampdata <- data.frame(subtype=unique(anno_col$Subtype))
p_alldata <- data.frame()
#allampdata <- data.frame()
for(g in drive_f){
   countdata <- data.frame()
   f <- function(x) sum( x >= 1)
   for(i in unique(anno_col$Subtype)){
       cluster <- rownames(subset(anno_col,Subtype == i))
       mut_f_sub <- mut_f[g,cluster]
       amp_count <- data.frame(amp=apply(mut_f_sub,1,f))
       amp_count$without_amp <- length(cluster) - amp_count$amp
       names(amp_count) <- c(paste(g,"amp",sep="_"),paste(g,"without_amp",sep="_"))
       rownames(amp_count) <- i
       countdata <- rbind(amp_count,countdata)
   }
   result <- fisher.test(countdata,workspace = 2e8)
  # print(paste(g,result$p.value,sep="_"))
   pdata <- data.frame(pvalue_amp = result$p.value)
   rownames(pdata) <- g  
   #print(pdata)

   allampdata <- cbind(allampdata,countdata) 
   p_alldata <- rbind(pdata,p_alldata)
}
write.table(t(allampdata[,2:ncol(allampdata)]),file="gene_amp_count.txt",sep="\t",row.names=T,quote=F)
write.table(p_alldata,file="gene_amp_count_pvalue.txt",sep="\t",row.names=T,quote=F)
p_alldata_sig <- subset(p_alldata,pvalue_amp < 0.05)
sig_amp_p <- data.frame()
for(s in rownames(p_alldata_sig)){
    data_a <- data.frame(matrix(0,nrow=2,ncol=2))
    rownames(data_a) <- c("with_OX","other")
    colnames(data_a) <- grep(paste(s,"_",sep=""),colnames(allampdata),value=T)
    allampdata_OX <- allampdata[grep("OXPHOS",rownames(allampdata),value=T),colnames(data_a)]
    allampdata_without_OX <- allampdata[setdiff(rownames(allampdata),grep("OXPHOS",rownames(allampdata),value=T)),colnames(data_a)]
    OX_count <- data.frame(count=colSums(allampdata_OX))
    without_OXcount <- data.frame(count=colSums(allampdata_without_OX))
    data_a[2,1] <- without_OXcount$count[1]
    data_a[2,2] <- without_OXcount$count[2]
    data_a[1,1] <- OX_count$count[1]
    data_a[1,2] <- OX_count$count[2]
    pvalue <- fisher.test(data_a,workspace = 2e8)
    sig_amp_pdata <- data.frame(sig_amp_OX_vs_others_pvalue=pvalue$p.value)
    rownames(sig_amp_pdata) <- s
    if(s == "PTEN"){
      print(data_a)
     print(sig_amp_pdata)
    }
    sig_amp_p <- rbind(sig_amp_pdata,sig_amp_p)
}
write.table(sig_amp_p,file="amp_sig_OX_vs_others_pvalue.txt",sep="\t",row.names=T,quote=F)
#loss p value
alllossdata <- data.frame(subtype=unique(anno_col$Subtype))
p_alldata_1 <- data.frame()
for(g in drive_f){
   countdata <- data.frame()
   f <- function(x) sum( x <= -1)
   for(i in unique(anno_col$Subtype)){
       cluster <- rownames(subset(anno_col,Subtype == i))
       mut_f_sub <- mut_f[g,cluster]
       loss_count <- data.frame(loss=apply(mut_f_sub,1,f))
       loss_count$without_loss <- length(cluster) - loss_count$loss
       names(loss_count) <- c(paste(g,"loss",sep="_"),paste(g,"without_loss",sep="_"))
       rownames(loss_count) <- i
       countdata <- rbind(loss_count,countdata)
   }
   result <- fisher.test(countdata,workspace = 2e8)
   #print(paste(g,result$p.value,sep="_"))
   pdata <- data.frame(pvalue_loss = result$p.value)
   rownames(pdata) <- g
   #print(pdata)

   alllossdata <- cbind(alllossdata,countdata)
   p_alldata_1 <- rbind(pdata,p_alldata_1)
}
write.table(t(alllossdata[,2:ncol(alllossdata)]),file="gene_loss_count.txt",sep="\t",row.names=T,quote=F)
write.table(p_alldata_1,file="gene_loss_count_pvalue.txt",sep="\t",row.names=T,quote=F)
p_alldata_1_sig <- subset(p_alldata_1,pvalue_loss < 0.05)
sig_loss_p <- data.frame()
for(s in rownames(p_alldata_1_sig)){   
     data_a <- data.frame(matrix(0,nrow=2,ncol=2))
     rownames(data_a) <- c("with_OX","other")
     colnames(data_a) <- grep(paste(s,"_",sep=""),colnames(alllossdata),value=T)
     alllossdata_OX <- alllossdata[grep("OXPHOS",rownames(alllossdata),value=T),colnames(data_a)]
     alllossdata_without_OX <- alllossdata[setdiff(rownames(alllossdata),grep("OXPHOS",rownames(alllossdata),value=T)),colnames(data_a)]
     OX_count <- data.frame(count=colSums(alllossdata_OX))
     without_OXcount <- data.frame(count=colSums(alllossdata_without_OX))
     data_a[2,1] <- without_OXcount$count[1]
     data_a[2,2] <- without_OXcount$count[2]
     data_a[1,1] <- OX_count$count[1]
     data_a[1,2] <- OX_count$count[2]
     pvalue <- fisher.test(data_a,workspace = 2e8)
     sig_loss_pdata <- data.frame(sig_loss_OX_vs_others_pvalue=pvalue$p.value)
    rownames(sig_loss_pdata) <- s
    sig_loss_p <- rbind(sig_loss_pdata,sig_loss_p)
}
write.table(sig_loss_p,file="loss_sig_OX_vs_others_pvalue.txt",sep="\t",row.names=T,quote=F)
p_gain_loss <- p_alldata
p_gain_loss$pvalue_loss <- p_alldata_1$pvalue_loss
display_numbers = matrix(ifelse(p_gain_loss >= 0.05, "",
    ifelse(p_gain_loss < 0.05 & p_gain_loss >= 0.01, "*",
        ifelse(p_gain_loss < 0.01 & p_gain_loss > 0.001, "**", "***"))), nrow(p_gain_loss))
rownames(display_numbers) <- rownames(p_gain_loss)
colnames(display_numbers) <- colnames(p_gain_loss)
display_numbers <- display_numbers[rownames(pat_mut_1),]
write.table(display_numbers,file="gene_amp_loss_pvalue_sig.txt",sep="\t",row.names=T,quote=F)
write.table(p_gain_loss,file="gene_amp_loss_pvalue.txt",sep="\t",row.names=T,quote=F)
pdf("cna_f.pdf",5,10)
p <- pheatmap(pat_mut_1,
              #scale="row",
              #color = colorRampPalette(colors = c("#87CEFA","#FA8072"))(2),legend_breaks=c(0,1),
              color = colorRampPalette(colors = c("navy", "white", "firebrick3"))(50),
              #legend=F,
              annotation_col = anno_col,
              #annotation_colors = ann_color,
              #annotation_row = annotation_row,
              cluster_cols=F,cluster_row=T,
              show_rownames = T, show_colnames = F)
print(p)
dev.off()
pdf("amp_sig_OX_vs_others_sig.pdf",5,10)
p <- pheatmap(pat_mut_1[rownames(subset(sig_amp_p,sig_amp_OX_vs_others_pvalue < 0.05)),],
              #scale="row",
              #color = colorRampPalette(colors = c("#87CEFA","#FA8072"))(2),legend_breaks=c(0,1),
              color = colorRampPalette(colors = c("navy", "white", "firebrick3"))(50),
              #legend=F,
              annotation_col = anno_col,
              main="OX_vs_others_amp_sig",
              #annotation_colors = ann_color,
              #annotation_row = annotation_row,
              cluster_cols=F,cluster_row=T,
              show_rownames = T, show_colnames = F)
print(p)
dev.off()
pdf("loss_sig_OX_vs_others_sig.pdf",5,12)
p <- pheatmap(pat_mut_1[rownames(subset(sig_loss_p,sig_loss_OX_vs_others_pvalue < 0.05)),],
              #scale="row",
              #color = colorRampPalette(colors = c("#87CEFA","#FA8072"))(2),legend_breaks=c(0,1),
              color = colorRampPalette(colors = c("navy", "white", "firebrick3"))(50),
              #legend=F,
              annotation_col = anno_col,
              main="OX_vs_others_loss_sig",
              #annotation_colors = ann_color,
              #annotation_row = annotation_row,
              cluster_cols=F,cluster_row=T,
              show_rownames = T, show_colnames = F)
print(p)
dev.off()
OX_vs_others_sig <- c(rownames(subset(sig_loss_p,sig_loss_OX_vs_others_pvalue < 0.05)),rownames(subset(sig_amp_p,sig_amp_OX_vs_others_pvalue < 0.05)))
length(OX_vs_others_sig)
length(unique(OX_vs_others_sig))
pdf("tnbc_OX_vs_others_amp_loss_sig.pdf",5,15)
p <- pheatmap(pat_mut_1[unique(OX_vs_others_sig),],
              #scale="row",
              #color = colorRampPalette(colors = c("#87CEFA","#FA8072"))(2),legend_breaks=c(0,1),
              color = colorRampPalette(colors = c("navy", "white", "firebrick3"))(50),
              #legend=F,
              annotation_col = anno_col,
              main="OX_vs_others_amp_loss_sig",
              #annotation_colors = ann_color,
              #annotation_row = annotation_row,
              cluster_cols=F,cluster_row=T,
              show_rownames = T, show_colnames = F)
print(p)
dev.off()


#.....................................TF activity in 4 subtypes...........................

library(viper)
library(magrittr)
options(stringsAsFactors=F)
setwd('./')
pfx <- 'LuminalB-TNBC_tumor'
outdir <- 'output_LuminalB-TNBC_tumor_ssviper/'
dir.create(outdir)

# create ExpressionSet
message(Sys.time(), ' :: create ExpressionSet... ')

# prepare input data --------------
## expression data QC
tnbc_all <- readRDS('./TPM.Rds')
str(tnbc_all)

rownames(tnbc_all) <- tnbc_all$GeneSymbol
tnbc_all$GeneSymbol <- NULL
str(tnbc_all)

### expr_gene_num > 5000 TNBC cel clusters
gene_num <- apply(tnbc_all, 2, function(x) { length(x[x>0]) })
str(gene_num)

gt_5k_sam <- gene_num[gene_num > 5000] %>% names()
str(gt_5k_sam)


### subset samples based on tumor purity
meta_data <- read.delim('./metadata_for-RNA.tsv',
    header = T, check.names = F) %>%
	subset((SampleName_RNA %in% gt_5k_sam) & (`Tumor (%)` > 75 | `normal epithelial` >= 90))
str(meta_data)

TN_samples <- subset(meta_data, `Tumor (%)` > 75 | `normal epithelial` >= 90) %>% .$SampleName_RNA
expr_TN <- tnbc_all[, TN_samples] %>% as.matrix
str(expr_TN)

### gene_qc
num_cell <- apply(expr_TN, 1, function(x) { length(x[x > 0]) })
ge_3_genes <- num_cell[num_cell >= 3] %>% names()
str(ge_3_genes)

expr_TN <- expr_TN[ge_3_genes, ]
expr_TN <- log2(expr_TN + 1)
str(expr_TN)

### celltype info
meta <- subset(meta_data, SampleName_RNA %in% colnames(expr_TN))
rownames(meta) <- meta$SampleName_RNA
meta$Celltype <- 'Tumour'
meta[meta$`Tumor (%)` < 75, 'Celltype'] <- 'Normal epithelial'
meta <- meta[, c('Celltype', 'Site', 'Patient')]
str(meta)


## build data
metadata <- data.frame(labelDescription=c('Celltype', 'Site', 'Patient'),
                       row.names=c('Celltype', 'Site', 'Patient'))
phenoData <- new("AnnotatedDataFrame", data = meta, varMetadata = metadata)
eset <- ExpressionSet(assayData = expr_TN, phenoData = phenoData)

# aracne2regulon
regulons <- readRDS('input/BRCA_network_regulon.Rds')

# prepare input for msviper
## signature
signature <- viperSignature(eset, "Celltype", "Normal epithelial", per = 1000,  verbose=T, seed = 1, cores = 10)
vpres <- viper(signature,  regulons, verbose = T, cores = 10)

save(eset, signature, vpres, file =  paste0(outdir, pfx, '_ssVIPER_result.Robj'))

# plot
d1 <- exprs(vpres)
colnames(d1) <- pData(vpres)$Site
dd <- dist(t(d1), method = "euclidean")

pdf(paste0(outdir, pfx, '_ssVIPER_similarity.pdf'), 20, 20)
heatmap(as.matrix(dd), Rowv = as.dendrogram(hclust(dd, method = "ward.D2")), symm = T)
dev.off()


#.....................................tumor subtype in single cell TNBC cohort...........................
library(Seurat)
library(GSVA)
library(dplyr)
library(survival)
library(survminer)
library(pheatmap)

info <- read.csv("1872-BIOKEY_metaData_cohort1_web.csv",header=T,row.names=1,stringsAsFactors=F)
tnbc_info <- subset(info,BC_type=="TNBC" & cellType=="Cancer_cell")

sc <- readRDS("pd1_add_meta_data.rds")
exp <- GetAssayData(object = sc, slot = "data")
matchcellname <- intersect(rownames(tnbc_info),colnames(exp))
tnbc_exp_match <- exp[,matchcellname]

tnbc_info_match <- tnbc_info[matchcellname,]

data <- t(as.matrix(tnbc_exp_match))
cellname<-rownames(data)
genename<-colnames(data)

ndata <- data[apply(data,1,function(x) !all(x==0)),]
tdata <- t(ndata)
tdata <- tdata[apply(tdata,1,function(x) !all(x==0)),]
tdata[1:3,1:3]

genelist <- readLines("./sub_genelist_gmt.txt")

flie=list()
for(i in 1:length(genelist)){
  genel <- strsplit( genelist[i],"\t")
  genechar <- as.vector(genel)
  gset<-c()
  for(ind in 2:length(genechar[[1]]) ){
    gset[ind-1] = genechar[[1]][(ind)]
  }
  listname=genechar[[1]][1]
  flie[[ genechar[[1]][1] ]] = gset
}

gbm_es <- gsva(tdata, flie, mx.diff=TRUE, verbose=FALSE)
#rownames(gbm_es) <- "OXPHOS"
t_gbm_es <- t(gbm_es)
write.table(t_gbm_es, file = "pd1_cohort1_subtype_gsva.tsv", sep = "\t",,quote=F)
anno_row <- tnbc_info_match[,c('timepoint','expansion','patient_id')]
pdf("pd1_cohort1_subtype_gsva.pdf",6,10)
p <- pheatmap(t_gbm_es,
              scale="row",cluster_col=F,
              annotation_row = anno_row,
              #cutree_rows=cutree_k,
              main="pd1_cohort1_subtype_gsva",
              show_colnames = T,show_rownames = F,clustering_method="ward.D2",
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               breaks = seq(-2,2,length.out = 50))
print(p)
dev.off()
#gsva_max <- data.frame(group=apply(t_gbm_es, 1, function(row) names(t_gbm_es)[which.max(row)]))

library(pheatmap)
library(dendextend)
library(magrittr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
info <- read.csv("1872-BIOKEY_metaData_cohort1_web.csv",header=T,row.names=1,stringsAsFactors=F)
tnbc_immu <- subset(info,BC_type=="TNBC" & cellType !="Cancer_cell")
tnbc_info <- subset(info,BC_type=="TNBC" & cellType=="Cancer_cell")

sc <- readRDS("pd1_add_meta_data.rds")
exp <- GetAssayData(object = sc, slot = "data")
matchcellname <- intersect(rownames(tnbc_info),colnames(exp))
tnbc_exp_match <- exp[,matchcellname]

tnbc_info_match <- tnbc_info[matchcellname,]
table(tnbc_info_match$cellType)

t_gbm_es <- read.table("pd1_cohort1_subtype_gsva.tsv",header=T,row.names=1)
t_gbm_es[1:5,]
identical(rownames(t_gbm_es),rownames(tnbc_info_match))

#anno_row <- tnbc_info_match[,c('timepoint','expansion','patient_id')]
anno_row <- tnbc_info_match[,c('timepoint','expansion')]

pdf("pd1_cohort1_subtype_gsva.pdf",4,5)
p <- pheatmap(t_gbm_es,
              scale="row",cluster_col=F,
              #cutree_rows=cutree_k,
              annotation_row = anno_row,
              main="pd1_cohort1_subtype_gsva",
              show_colnames = T,show_rownames = F,clustering_method="ward.D2",
              #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              color = colorRampPalette(c('#137AB5','white','#CC0724'))(50),
               breaks = seq(-2,2,length.out = 50))
print(p)
dev.off()

cutree_k <- 16 
pdf("pd1_cohort1_subtype_gsva_cutree.pdf",5,5)
p <- pheatmap(t_gbm_es,
              scale="row",cluster_col=F,
              cutree_rows=cutree_k,
              annotation_row = anno_row,
              main="pd1_cohort1_subtype_gsva_cutree",
              show_colnames = T,show_rownames = F,clustering_method="ward.D2",
              #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              color = colorRampPalette(c('#137AB5','white','#CC0724'))(50),
               breaks = seq(-2,2,length.out = 50))
print(p)
dev.off()

dend <- p$tree_row %>%
        as.dendrogram() %>%
        set("labels_to_character") %>%
        color_branches(k = cutree_k, groupLabels = TRUE)

pdf("pd1_cohort1_subtype_gsva_cutree_dend.pdf", 4, 10)
plot(rev(dend), main = "Row dendrogram", horiz=T)
dev.off()

clusters <- dendextend::cutree(dend, cutree_k, order_clusters_as_data = FALSE)
clusters_ex <- data.frame(tree=dendextend::cutree(dend, cutree_k, order_clusters_as_data = F))
clusters_ex$group <- "subtype"
clusters_ex[which(clusters_ex$tree == "1"),]$group <- "OXPHOS"
clusters_ex[which(clusters_ex$tree == "2"),]$group <- "OXPHOS+IM/MES"
clusters_ex[which(clusters_ex$tree == "3"),]$group <- "OXPHOS"
clusters_ex[which(clusters_ex$tree == "4"),]$group <- "OXPHOS"
clusters_ex[which(clusters_ex$tree == "5"),]$group <- "IM regulatory+OXPHOS"
clusters_ex[which(clusters_ex$tree == "6"),]$group <- "IM regulatory"
clusters_ex[which(clusters_ex$tree == "7"),]$group <- "IM regulatory"
clusters_ex[which(clusters_ex$tree == "8"),]$group <- "LAR"
clusters_ex[which(clusters_ex$tree == "9"),]$group <- "IM regulatory"
clusters_ex[which(clusters_ex$tree == "10"),]$group <- "IM regulatory"
clusters_ex[which(clusters_ex$tree == "11"),]$group <- "IM regulatory"
clusters_ex[which(clusters_ex$tree == "12"),]$group <- "IM/MES"
clusters_ex[which(clusters_ex$tree == "13"),]$group <- "IM regulatory"
clusters_ex[which(clusters_ex$tree == "14"),]$group <- "IM regulatory"
clusters_ex[which(clusters_ex$tree == "15"),]$group <- "IM regulatory+OXPHOS"
clusters_ex[which(clusters_ex$tree == "16"),]$group <- "IM regulatory"

head(clusters_ex)
tnbc_info <- tnbc_info_match[rownames(clusters_ex),]

clusters_ex[,3:6] <- tnbc_info[,c(3:5,8)]
clusters_ex$cellbarcode <- rownames(clusters_ex)
clusters_ex$group <- gsub("IM regulatory","IM",clusters_ex$group)
clusters_ex$patient <- paste(clusters_ex$patient_id,clusters_ex$expansion,sep="_")

a <- data.frame(table(clusters_ex$patient,clusters_ex$group))
a <- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100)
pat_group_percent <- dcast(a,Var1~Var2,value.var = "percent")
rownames(pat_group_percent) <- pat_group_percent$Var1
pat_group_percent$Var1 <- NULL
pat_group_percent <- pat_group_percent[rev(rownames(pat_group_percent)),rev(colnames(pat_group_percent))]
write.table(pat_group_percent,file="pd1_cohort1_pat_group_percent.tsv",sep="\t",row.names=T,quote=F)

tnbc_immu$patient <- paste(tnbc_immu$patient_id,tnbc_immu$expansion,sep="_")

b <- data.frame(table(tnbc_immu$patient,tnbc_immu$cellType))
b <- ddply(b,.(Var1),transform,percent=Freq/sum(Freq)*100)
pat_celltype_percent <- dcast(b,Var1~Var2,value.var = "percent")
rownames(pat_celltype_percent) <- pat_celltype_percent$Var1
pat_celltype_percent$Var1 <- NULL
pat_celltype_percent <- pat_celltype_percent[rev(rownames(pat_celltype_percent)),rev(colnames(pat_celltype_percent))]
write.table(pat_celltype_percent,file="pd1_cohort1_pat_celltype_percent.tsv",sep="\t",row.names=T,quote=F)

#data <- clusters_ex  %>%
 # group_by(group) %>%
 # mutate(Total = nrow(clusters_ex)) %>%
 # ungroup()

pdf("pd1_cohort1_subgroup_percent.pdf",13,30)
p0 <- ggplot(clusters_ex,aes(x=cohort,fill=group))+geom_bar(position="fill")+ylab("Cell (%)")+theme_classic()+
      #scale_fill_manual(values=colors)+
      scale_y_continuous(expand = c(0,0),labels = c(0,25,50,75,100))+
      ggtitle("Group_percent")+
      theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+ theme_bw(base_size=20)

p1 <- ggplot(clusters_ex,aes(x=group,fill=patient_id))+geom_bar(position="fill")+ylab("Cell (%)")+theme_classic()+
      #scale_fill_manual(values=colors)+
      scale_y_continuous(expand = c(0,0),labels = c(0,25,50,75,100))+
      ggtitle("Group_patient_percent")+
      theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+ theme_bw(base_size=20)
      #+labs(fill = "Invasiveness")

aa <- as.data.frame.matrix(table(clusters_ex$group,clusters_ex$timepoint))
row_sums <- rowSums(aa)
df_proportions <- sweep(aa, 1, row_sums, FUN = "/")
chi_square_test <- chisq.test(df_proportions)
pvalue1 <- round(chi_square_test$p.value,3)

p2 <- ggplot(clusters_ex,aes(x=group,fill=timepoint))+geom_bar(position="fill")+ylab("Cell (%)")+theme_classic()+
      #scale_fill_manual(values=colors)+
      scale_y_continuous(expand = c(0,0),labels = c(0,25,50,75,100))+
      ggtitle(paste("Group_timepoint_percent","pvalue=",pvalue1,sep="  "))+
      theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+ theme_bw(base_size=20)

clusters_ex_E <- subset(clusters_ex,expansion=="E")
dim(clusters_ex_E)
bb <- as.data.frame.matrix(table(clusters_ex_E$timepoint,clusters_ex_E$group))
row_sums_bb <- rowSums(bb)
row_sums_bb
df_proportions_bb <- sweep(bb, 1, row_sums_bb, FUN = "/")
chi_square_test_bb <- chisq.test(df_proportions_bb)
pvalue2 <- round(chi_square_test_bb$p.value,3)

p3 <- ggplot(clusters_ex_E,aes(x=timepoint,fill=group))+geom_bar(position="fill")+ylab("Cell (%)")+theme_classic()+
      #scale_fill_manual(values=colors)+
      scale_y_continuous(expand = c(0,0),labels = c(0,25,50,75,100))+
      ggtitle(paste("E_timepoint_group_percent","pvalue=",pvalue2,sep="  "))+
      theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+ theme_bw(base_size=20)

clusters_ex_NE <- subset(clusters_ex,expansion=="NE")
dim(clusters_ex_NE)
cc <- as.data.frame.matrix(table(clusters_ex_NE$timepoint,clusters_ex_NE$group))
row_sums_cc <- rowSums(cc)
df_proportions_cc <- sweep(cc, 1, row_sums_cc, FUN = "/")
df_proportions_cc
chi_square_test_cc <- chisq.test(df_proportions_cc)
pvalue3 <- round(chi_square_test_cc$p.value,3)

p4 <- ggplot(clusters_ex_NE,aes(x=timepoint,fill=group))+geom_bar(position="fill")+ylab("Cell (%)")+theme_classic()+
      #scale_fill_manual(values=colors)+
      scale_y_continuous(expand = c(0,0),labels = c(0,25,50,75,100))+
      ggtitle(paste("NE_timepoint_group_percent","pvalue=",pvalue3,sep="  "))+
      theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+ theme_bw(base_size=20)

clusters_ex$patient <- paste(clusters_ex$patient_id,clusters_ex$expansion,sep="_")
p5 <- ggplot(clusters_ex,aes(x=patient,fill=group))+geom_bar(position="fill")+ylab("Cell (%)")+theme_classic()+
      #scale_fill_manual(values=colors)+
      scale_y_continuous(expand = c(0,0),labels = c(0,25,50,75,100))+
      ggtitle("patient_group_percent")+
      theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+ theme_bw(base_size=20)+coord_flip()

tnbc_immu$patient <- paste(tnbc_immu$patient_id,tnbc_immu$expansion,sep="_")
p6 <- ggplot(tnbc_immu,aes(x=patient,fill=cellType))+geom_bar(position="fill")+ylab("Cell (%)")+theme_classic()+
      #scale_fill_manual(values=colors)+
      scale_y_continuous(expand = c(0,0),labels = c(0,25,50,75,100))+
      ggtitle("patient_cellType_percent")+
      theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+ theme_bw(base_size=20)+coord_flip()

combined_plot <- p0 / p1 / p2 / p3 / p4 / p5 / p6
print(combined_plot)
dev.off()


#.....................................correlation of cellular type in group E...........................
library(pheatmap)
getwd()
setwd('./')
df= read.delim("proportion-of-cell-type.txt",header = T,row.names = 1,check.names=FALSE)   
file<-read.table("cor-order.txt",header=T,row.names=1,check.names=FALSE)
desired_order<-file[,2]

# 检查列名是否匹配
if (!all(desired_order %in% colnames(df))) {
  stop("列名不匹配，请检查列名是否正确")
}
df<-df[,desired_order]
r <- cor(df,method = "spearman") 
r[lower.tri(r)] <- NA
pheatmap(r, 
         show_colnames = TRUE,   
         show_rownames=TRUE,     
         fontsize=5,             
         color=c(colorRampPalette(colors = c("#137AB5","white"))(50),colorRampPalette(colors = c("white","#CC0724"))(50)),
         annotation_legend=TRUE, 
         border_color=NA,        
         scale="none",           
         cluster_rows = FALSE,    
         cluster_cols = FALSE,     
         main="cohort1-Pre-NE-QC400-author-annotate-reorder",
         na_col="white"
)

