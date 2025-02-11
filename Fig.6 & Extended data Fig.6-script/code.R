#script for figure6
#.....................................gene expression clustering...........................
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(analogue)
library(pheatmap)
library(patchwork)
# input -----------
## expression
expr_all <- readRDS('./TPM.Rds') %>%
  data.frame(row.names = 1, check.names = F) %>% select(matches("A|D|BC"))
#%>%
  #{ log2(. + 1) }
#%>% select(matches("BC"))
num=ncol(expr_all)
print(paste('total cell cluster: ',num))
meta_data <- read.delim('./metadata.tsv',header = T, check.names = F)
#meta_data<-subset(meta_data,`Site`=='primary')
meta_data<-subset(meta_data,`Tumor (%)`>75)
meta_data$SampleName_purity <- paste0(meta_data$SampleName, '-',
                                      meta_data$`Tumor (%)`, '%')
colnames(expr_all) %<>% plyr::mapvalues(from = meta_data$SampleName_RNA,
                                        to = meta_data$SampleName_purity)
expr_all<-expr_all %>% select(contains("_"))
rb.gene <- rownames(expr_all)[grep('^RP[SL]',rownames(expr_all))]
expr_all <- expr_all[setdiff(rownames(expr_all),rb.gene),]
## QC
### expr_gene_num > 5000 TNBC cel clusters
gene_num <- apply(expr_all, 2, function(x) { length(x[x>0]) })
gt_5k_sam <- gene_num[gene_num > 5000] %>% names()
str(gt_5k_sam)
expr_all <- expr_all[, gt_5k_sam] %>% as.data.frame
print(head(expr_all))
print(paste('after QC, cell cluster: ',ncol(expr_all)))

### Find which genes are >=5TPM in >=5%  of tumour samples
perc5 <- ceiling(ncol(expr_all)*.05)
atl10 <- apply(expr_all,1,function(x) sum(x>=20))
#atl10 <- apply(expr_all,1,function(x) sum(x>=5))
expr_all$GeneID <- rownames(expr_all)
filens <- expr_all$GeneID[which(atl10>=perc5)]
expr <- expr_all[filens,]
expr_f <- expr[,1:(ncol(expr)-1)]
expr_f_log2 <- log2(expr_f + 1)
#PID <- unique(sub('_.*', '', colnames(expr_f)))
PID <- unique(sub('_[1-9].*', '', colnames(expr_f)))
print(PID)
## calculate mean_expression and sd for each gene
mean_exp <- data.frame()
sd_exp <- data.frame()
for(i in PID){
   pid <- grep(i,colnames(expr_f),value=T)
   print(i)
  # print(pid)
   exp_pid <- expr_f[,pid]
   PID_exp_mean <- data.frame(meanexp = log2(rowMeans(exp_pid) + 1))
   #PID_exp_sd <- data.frame(sdexp = rowSds(exp_pid))
   PID_exp_sd <- data.frame(sdexp = apply(expr_f_log2[, pid],1,sd))
   colnames(PID_exp_mean) <- i
   colnames(PID_exp_sd) <- i
   if(i=='P02_BC'){
      mean_exp <- PID_exp_mean
      sd_exp <- PID_exp_sd
   }else{
      mean_exp <- cbind(mean_exp,PID_exp_mean)
      sd_exp <- cbind(sd_exp,PID_exp_sd)
   }
}
# update patient id
id <- read.delim('./tumor_id.txt')
colnames(mean_exp) %<>% plyr::mapvalues(from = id$Origin_ID, to = id$New_ID) 
colnames(sd_exp)  %<>% plyr::mapvalues(from = id$Origin_ID, to = id$New_ID)
# cluster genes ----------
dexp <- dist(mean_exp, method='euclidean')
dsd <- dist(sd_exp[rownames(mean_exp), colnames(mean_exp)], method='euclidean')
dd <- fuse(dexp, dsd, weights=c(0.5,0.5))
clustord <- hclust(dd, method='ward.D2')
## order samples
dsamexp <- dist(t(mean_exp), method='euclidean')
dsamsd <- dist(t(sd_exp[rownames(mean_exp), colnames(mean_exp)]), method='euclidean')
ddsam <- fuse(dsamexp,dsamsd,weights=c(0.5,0.5))
clustsam <- hclust(ddsam,method='ward.D2')
samord <- colnames(sd_exp)[clustsam$order]
#mean_exp_reorder <- mean_exp[geneord, samord]
#sd_exp_reorder <- sd_exp[geneord, samord]
## set cluster id based on median of mean_mean_exp ( of genes included)
genedf <- data.frame(row.names = rownames(mean_exp),
	Mean_mean_exp = rowMeans(mean_exp),
	Mean_sd_exp = rowMeans(sd_exp)
)
k <- 4
tmpmemb <- cutree(clustord,k = k)
medstat <- boxplot(genedf$Mean_mean_exp~tmpmemb,plot=F)$stats[3,]
ordgroup <- order(medstat,decreasing=T)
genedf$Group <- rep('1',nrow(genedf))
genedf[tmpmemb==ordgroup[2],'Group'] <- '2'
genedf[tmpmemb==ordgroup[3],'Group'] <- '3'
genedf[tmpmemb==ordgroup[4],'Group'] <- '4'
## order matrix for visualization
geneord <- genedf[order(genedf$Group, -genedf$Mean_sd_exp), ] %>% row.names
mean_exp_reorder <- mean_exp[geneord, samord]
sd_exp_reorder <- sd_exp[geneord, samord]
# visual --------------------
p1 <- pheatmap(mean_exp_reorder,
	scale = 'column',
	cluster_rows=F,
	cluster_cols=F,
	show_rownames=F,
	main='\nMean expression',
	border = F,
	annotation_row = genedf,
        color = c(colorRampPalette(colors = c("#137AB5","white"))(50),colorRampPalette(colors = c("white","#CC0724"))(50)),
	#color = colorRampPalette(rev(brewer.pal(9, 'Spectral')))(100),
	breaks = unique(c(seq(-2,2, length=100))),
	cellwidth=15)
p2 <- pheatmap(sd_exp_reorder, 
	scale = 'column',
	cluster_rows=F,
	cluster_cols=F,
	show_rownames=F,
	main='\nStandard deviation of expression',
	annotation_row = genedf,
	border = F,
        color = c(colorRampPalette(colors = c("#137AB5","white"))(50),colorRampPalette(colors = c("white","#CC0724"))(50)),
	#color = colorRampPalette(rev(brewer.pal(9, 'Spectral')))(100),
	breaks = unique(c(seq(-3,3, length=100))),
	cellwidth=15)
dev.off()
pdf(paste0('./00.gene_clustering_heatmap_', k, 'group.pdf'), 15)
cowplot::plot_grid(p1$gtable,p2$gtable,ncol=2)
dev.off()
save(genedf, mean_exp_reorder, sd_exp_reorder, file = paste0('./00.gene_clustering_', k, 'group.Robj'))
write.csv(genedf,file=paste0('./00.gene_clustering_',".csv"))


#.....................................pathway expression clustering...........................
library(GSVA)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(analogue)
library(pheatmap)
library(patchwork)
#library(msigdbr)
# input -----------
## expression
expr_all <- readRDS('./TPM.Rds') %>%
  data.frame(row.names = 1, check.names = F) %>% select(matches("A|D"))
  #{ log2(. + 1) }
meta_data <- read.delim('./metadata.tsv',
                        header = T, check.names = F)
meta_data<-subset(meta_data,`Tumor (%)`>75)
meta_data$SampleName_purity <- paste0(meta_data$SampleName, '-',
                                      meta_data$`Tumor (%)`, '%')

colnames(expr_all) %<>% plyr::mapvalues(from = meta_data$SampleName_RNA,
                                        to = meta_data$SampleName_purity)
expr_all<-expr_all %>% select(contains("_"))
rb.gene <- rownames(expr_all)[grep('^RP[SL]',rownames(expr_all))]
expr_all <- expr_all[setdiff(rownames(expr_all),rb.gene),]
## QC
### expr_gene_num > 5000 TNBC cel clusters
gene_num <- apply(expr_all, 2, function(x) { length(x[x>0]) })
gt_5k_sam <- gene_num[gene_num > 5000] %>% names()
str(gt_5k_sam)
expr_all <- expr_all[, gt_5k_sam] %>% as.data.frame
print(head(expr_all))
print(paste('after QC, cell cluster: ',ncol(expr_all)))
## QC
### Find which genes are >=5TPM in >=5%  of tumour samples
perc5 <- ceiling(ncol(expr_all)*.05)
atl10 <- apply(expr_all,1,function(x) sum(x>=20))
#atl10 <- apply(expr_all,1,function(x) sum(x>=5))
expr_all$GeneID <- rownames(expr_all)
filens <- expr_all$GeneID[which(atl10>=perc5)]

expr <- expr_all[filens,]
expr_f <- expr[,1:(ncol(expr)-1)]
expr_f_log2 <- log2(expr_f + 1)

PID <- unique(sub('_[1-9].*', '', colnames(expr_f)))
print(PID)
# hallmark pathway
#m_df = msigdbr(species = "Homo sapiens", category = "H")
#Hallmark_cancer_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
Hallmark_cancer_list <- readRDS('./Hallmark_cancer_list.Rds')
names(Hallmark_cancer_list) <- gsub('HALLMARK_(\\S+)','\\1',names(Hallmark_cancer_list))
# Remove spermatogenesis, myogenesis and pancreas_beta_cells
Hallmark_cancer_list$SPERMATOGENESIS <- NULL
Hallmark_cancer_list$MYOGENESIS <- NULL
Hallmark_cancer_list$PANCREAS_BETA_CELLS <- NULL
# Rename COMPLEMENT
Hallmark_cancer_list$COMPLEMENT_INNATE_IMMUNE_SYSTEM <- Hallmark_cancer_list$COMPLEMENT
Hallmark_cancer_list$COMPLEMENT <- NULL
# Include WNT_signalling
geneinfo <- data.table::fread('./complete_gene_info.txt.gz')
Hallmark_cancer_list$WNT_SIGNALING <- read.table('./wnt_signalling_entrez.txt')[,1] %>%
  intersect(., geneinfo$Entrez) %>%
  plyr::mapvalues(from = geneinfo$Entrez, to = geneinfo$Symbol)
#unique gene in hallmarkgeneset
Hallmark_cancer_list<- lapply(Hallmark_cancer_list, unique)
## Run ssGSEA and calculate mean enrichment and mean standard deviation per pathway
expmat <- as.data.frame(matrix(0L,nrow=length(names(Hallmark_cancer_list)),ncol=length(PID)))
row.names(expmat) <- names(Hallmark_cancer_list)
colnames(expmat) <- PID
sdmat <- expmat
esl <- list()
for(i in PID) {
  print(i)
  pid <- grep(i, colnames(expr_f), value=T)
  exp_pid <- expr_f_log2[, pid]
  hallmark_gsva <- gsva(data.matrix(exp_pid),method='ssgsea',
                        Hallmark_cancer_list,verbose=FALSE)
  expmat[, i] <- rowMeans(as.matrix(hallmark_gsva))
  sdmat[, i] <- apply(hallmark_gsva,1,function(x) sd(x))
  esl[[i]] <- hallmark_gsva
}
es_df <- do.call(cbind, esl)
write.csv(es_df,file='./es_df.csv')
pathdf <- data.frame(Mean_Var=rowMeans(as.matrix(sdmat)),Mean_Mean_Enrichment=rowMeans(as.matrix(log(expmat+1))))
row.names(pathdf) <- names(Hallmark_cancer_list)
pathdf <- pathdf[which(pathdf$Mean_Mean_Enrichment>0),]
## Calculate distance matrices and cluster
# First filter for useable pathways
mean_exp <- expmat[row.names(pathdf),]
sd_exp <- sdmat[row.names(pathdf),]
# cluster genes ----------
dexp <- dist(mean_exp, method='euclidean')
dsd <- dist(sd_exp[rownames(mean_exp), colnames(mean_exp)], method='euclidean')
dd <- fuse(dexp, dsd, weights=c(0.5,0.5))
clustord <- hclust(dd, method='ward.D2')
## order samples
dsamexp <- dist(t(mean_exp), method='euclidean')
dsamsd <- dist(t(sd_exp[rownames(mean_exp), colnames(mean_exp)]), method='euclidean')
ddsam <- fuse(dsamexp,dsamsd,weights=c(0.5,0.5))
clustsam <- hclust(ddsam,method='ward.D2')
samord <- colnames(sd_exp)[clustsam$order]
## set cluster id based on median of mean_mean_exp ( of genes included)
genedf <- data.frame(row.names = rownames(mean_exp),
	Mean_mean_exp = rowMeans(mean_exp),
	Mean_sd_exp = rowMeans(sd_exp)
)
k <- 4
tmpmemb <- cutree(clustord,k = k)
medstat <- boxplot(genedf$Mean_mean_exp~tmpmemb,plot=F)$stats[3,]
ordgroup <- order(medstat,decreasing=T)
genedf$Group <- rep('1',nrow(genedf))
genedf[tmpmemb==ordgroup[2],'Group'] <- '2'
genedf[tmpmemb==ordgroup[3],'Group'] <- '3'
genedf[tmpmemb==ordgroup[4],'Group'] <- '4'
## order matrix for visualization
geneord <- genedf[order(genedf$Group, -genedf$Mean_sd_exp), ] %>% row.names
mean_exp_reorder <- mean_exp[geneord, samord]
sd_exp_reorder <- sd_exp[geneord, samord]
# visual --------------------
p1 <- pheatmap(mean_exp_reorder,
	scale = 'column',
	cluster_rows=F,
	cluster_cols=F,
	show_rownames=F,
	main='\nMean expression',
	border = F,
	annotation_row = genedf,
	color = c(colorRampPalette(colors = c("#137AB5","white"))(50),colorRampPalette(colors = c("white","#CC0724"))(50)),
	#color = colorRampPalette(rev(brewer.pal(9, 'Spectral')))(100),
	breaks = unique(c(seq(-2,2, length=100))),
	cellwidth=15)
p2 <- pheatmap(sd_exp_reorder,
	scale = 'column',
	cluster_rows=F,
	cluster_cols=F,
	show_rownames=T,
	main='\nStandard deviation of expression',
	annotation_row = genedf,
	border = F,
	color = c(colorRampPalette(colors = c("#137AB5","white"))(50),colorRampPalette(colors = c("white","#CC0724"))(50)),
	#color = colorRampPalette(rev(brewer.pal(9, 'Spectral')))(100),
	breaks = unique(c(seq(-3,3, length=100))),
	cellwidth=15)
dev.off()
pdf(paste0('./00.pathway_clustering_heatmap_', k, 'group.pdf'), 20)
cowplot::plot_grid(p1$gtable,p2$gtable,ncol=2)
dev.off()
save(genedf, es_df, mean_exp_reorder, sd_exp_reorder, file = paste0('./00.pathway_clustering_', k, 'group.Robj'))
write.csv(genedf,file=paste0('./00.gene_clustering_',".csv"))


#.....................................phylogenetic signals...........................
library(phytools)
library(dplyr)
library(magrittr)
library(pheatmap)
# input ---------------
## expr matrix
load('./00.pathway_clustering_4group.Robj')
es_df <- as.data.frame(es_df)
## gene filtering
gene_sel <- subset(genedf, Group != '5') %>% rownames
str(gene_sel)
#gene_sel <- rownames(genedf)
# call phylogenetic signals ------------ 
nwkDir <- "./input/DNA_tree_75TP_hclust"
files <- list.files(nwkDir, pattern="*75%.nwk", full.names=T)
PID <- unique(sub('_.*', '', basename(files)))
for(i in PID){
   print(paste0('Run for ', i))
   nwkfile <- grep(i, files, value=T)
   tree <- read.newick(nwkfile)
   exp_pid_t <- as.data.frame(t(es_df[intersect(tree$tip.label, colnames(es_df))]))
   treetrim <- ape::drop.tip(tree, tip = setdiff(tree$tip.label, rownames(exp_pid_t)))

   sig_df <- data.frame(Gene = gene_sel, Lambda = NA, Lpval = NA)

   for(j in gene_sel) {

 	  expr_named <- setNames(exp_pid_t[, j], rownames(exp_pid_t))
	  
	  if (var(expr_named) > 0) {
		res <- phylosig(treetrim, expr_named, method="lambda", test=T)
        sig_df[sig_df$Gene == j,'Lambda'] <- res$lambda
		sig_df[sig_df$Gene == j,'Lpval'] <- res$P
      } else {
    	res <- phylosig(treetrim, expr_named, method='lambda', test=F)
        sig_df[sig_df$Gene == j,'Lambda'] <- res$lambda
		sig_df[sig_df$Gene == j,'Lpval'] <- 1
      }

   }    
   write.csv(sig_df, file=paste0('./01.', i, "_sigphylo_pathway_Pvalue.csv"), quote=F, row.names=F)
}


#.......................................meta-KEGG analysis..................................
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
input_path='./input'
out_path='./output'
files <- list.files(paste0(input_path,'/input/'),pattern=".txt")
#dir.create(out_path)
for(i in files){
   pid <- sub(".txt","",i)
   print(i)
   genes <- read.table(paste0(input_path,'/input/',i),header=T,stringsAsFactors=F)
   print('step1')
   #把symbol转成gene_id
   gene2ID = bitr(genes$x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
   ek <- enrichKEGG(gene = gene2ID$ENTREZID, #需要分析的基因的EntrezID
               organism = "hsa", #人基因数据库
               keyType = "kegg",
               #pvalueCutoff = 0.05, #设置pvalue界值
               pAdjustMethod = "fdr",
               qvalueCutoff = 0.05, #设置qvalue界值(FDR校正后的p值）
               use_internal_data = FALSE)
   #print('step1: convert')
   # ek2 = setReadable(ek, #前面分析的结果
   #                OrgDb = "org.Hs.eg.db", #人类数据库
   #                keyType = "ENTREZID") #要转换的基因类型
   write.table(ek,file=paste0(out_path,'/',pid,"_ek.txt"), sep="\t", quote=F, row.names = F)
   
   if(nrow(ek) != 0){
     pdf(file=paste0(out_path,'/',pid,"_kegg_barplot.pdf"),width = 12,height = 10)
     p <-  barplot(ek, x = "GeneRatio", color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
                   showCategory =15) #只显示前10
     print(p)
     dev.off()
     
     print('step3')
     pdf(file=paste0(out_path,'/',pid,"_kegg_dotplot.pdf"),width = 12,height = 10)
     p1 <- dotplot(ek,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
                   showCategory =15)#只显示前5#以ONTOLOGY类型分屏绘图
     print(p1)
     dev.off()
   }
   
   eG <- enrichGO(gene = gene2ID$ENTREZID, #需要分析的基因的EntrezID
                  OrgDb = org.Hs.eg.db, #人基因数据库
                  #pvalueCutoff =0.05, #设置pvalue界值
                  pAdjustMethod = "fdr",
                  qvalueCutoff = 0.05, #设置qvalue界值(FDR校正后的p值）
                  ont="all") #选择功能富集的类型，可选BP、MF、CC，这里选择all。
   write.table(eG,file=paste0(out_path,'/',pid,"_eG.txt"), sep="\t", quote=F, row.names = F)
   
   
   # print('step2')
   if(nrow(eG) != 0){
     pdf(file=paste0(out_path,'/',pid,"_go_barplot.pdf"),width = 12,height = 10)
     p <-  barplot(eG, x = "GeneRatio", color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
        showCategory =15,cex.axis=20,cex.names=20) #只显示前10
     print(p)
     dev.off()

     print('step3')
     pdf(file=paste0(out_path,'/',pid,"_go_dotplot.pdf"),width = 12,height = 10)
     p1 <- dotplot(eG,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =15)#只显示前5#以ONTOLOGY类型分屏绘图
     print(p1)
     dev.off()
   }
   
   print('step4')
   #enricher
   gene_name<-gene2ID$SYMBOL
   gene_ID<-gene2ID$ENTREZID
   annotation<-read.table(paste0('./','/enricher-input_level2.txt'),header=T,check.names=FALSE)
   term2gene=data.frame(Term=annotation$pathway_id, Gene=annotation$ko_ID)
   term2name=data.frame(Term=annotation$pathway_id,Name=annotation$pathway_name)
   ENRICH=enricher(gene=gene_ID,
                   pAdjustMethod='fdr',
                   qvalueCutoff = 0.1,
                   TERM2GENE=term2gene,
                   TERM2NAME=term2name
                   )
   write.table(ENRICH,file=paste0(out_path,'/',pid,"_enricher-level2.txt"), sep="\t", quote=F, row.names = F)
   gene_name<-gene2ID$SYMBOL
   gene_ID<-gene2ID$ENTREZID
   annotation<-read.table(paste0('./','/enricher-input_level3.txt'),header=T,check.names=FALSE)
   term2gene=data.frame(Term=annotation$pathway_id, Gene=annotation$ko_ID)
   term2name=data.frame(Term=annotation$pathway_id,Name=annotation$pathway_name)
   ENRICH=enricher(gene=gene_ID,
                   pAdjustMethod='fdr',
                   qvalueCutoff = 0.1,
                   TERM2GENE=term2gene,
                   TERM2NAME=term2name
   )
   write.table(ENRICH,file=paste0(out_path,'/',pid,"_enricher-level3.txt"), sep="\t", quote=F, row.names = F)
}


#.....................................cophenetic correlation.......................
library(pheatmap)
library(wesanderson)
library(dendextend)
library(magrittr)
library(ape)
options(datatable.fread.datatable=FALSE)
rm(list =  ls())
setwd('./')
#method <- 'ward.D2'
method <- 'average'
outdir <- paste0('TNBC_output_', method, '/')
dir.create(outdir)
regcol <- c(A='#E31A1C',BC='#377DB8',D='#4DAE49', Root='#808080')
rcol <- regcol[c(1:3)]
## genes for clustering
expr <- load("./00.pathway_clustering_4group.Robj")
es_df[1:5,1:5]
## rename samples
meta_data <- read.delim('./metadata.tsv',
                        header = T, check.names = F)
meta_data$SampleName_purity <- paste0(meta_data$SampleName, '-',
                                      meta_data$`Tumor (%)`, '%')
colnames(es_df) %<>% plyr::mapvalues(from = meta_data$SampleName_RNA, to = meta_data$SampleName_purity)
# tree data
df <- data.frame(Patient = sort(unique(sub('_.*', '', colnames(es_df)))), Cophenetic_corr = '-')
for (pat in setdiff(unique(sub('_.*', '', colnames(es_df))),c('ZLR'))) {
   gene_sel <- rownames(es_df)
   if (length(gene_sel) > 0) {
  tre <- read.tree(paste0('./DNA_tree_75TP_kmeans/', pat, '_tree_75%.nwk'))
  shared_sam <- intersect(colnames(es_df), tre$tip.label)
  trimtree <- drop.tip(tre, tip = setdiff(tre$tip.label, c(shared_sam))) #, 'Normal')))
  ## Plot DNA phylogeny
  pdf(paste0(outdir, pat, '_DNA-pathway_cophenetic-corr.pdf'))
  regions <- gsub('(.*)_(.*)_(.*)','\\2', trimtree$tip.label)
  options(scipen = -1);par(mar=c(0,0,2,0),xpd=T)
  plot.phylo(trimtree, type="phylogram", align.tip.label = T,edge.width=3,
             font=2,cex=0.9, tip.color=regcol[c(regions,'Root')], label.offset=20)
  title(main=pat)
  ## Remove gene with no variance (i.e. zero expression in particular tumour)
  gene_expr <- as.data.frame(t(es_df[gene_sel, shared_sam]))
  curvars <- apply(gene_expr, 2, sd)
  explot <- gene_expr[, which(curvars != 0)]
  # Plot clustered heatmap
  rowannot <- data.frame(Region=gsub('(.*)_(.*)_(.*)','\\2', row.names(explot)))
  row.names(rowannot) <- row.names(explot)
  explot_scaled <- scale(explot)
  explot_scaled[explot_scaled > 3] <- 3 
  explot_scaled[explot_scaled < -3] <- -3
  mybreaks <- seq(floor(min(explot_scaled)),ceiling(max(explot_scaled)),by=0.05)
  mat <- pheatmap(explot_scaled,show_rownames = T,show_colnames=F,cluster_rows=T,cluster_cols=T,treeheight_col=0,
                  border_color=NA,fontsize = 10,legend = T,treeheight_row = 100, scale = 'none', clustering_method = method,
                  color=wes_palette("Zissou1", length(mybreaks)-1, type = "continuous"), breaks=mybreaks)
  ## Plot DNA phylogenies as dendrograms, clustered heatmaps and get tanglegram statistics
  ### Extended Data Figure 5 (deconstructed)
  patdend <- as.dendrogram(mat$tree_row)
  patdend <- set(patdend, "labels_color", value = rcol[gsub('(.*)_(.*)_(.*)', '\\2', labels(patdend))])
  patdend <- set(patdend,"branches_lwd", 3)
  dnadend <- trimtree %>% phylogram::as.dendrogram.phylo()
  dnadend <- set(dnadend, "labels_color", value = rcol[gsub('(.*)_(.*)_(.*)', '\\2', labels(dnadend))])
  dnadend <- set(dnadend, "branches_lwd", 3)
  bothdends <- dendlist(dnadend, patdend)
  curcor <- cor.dendlist(bothdends)
  curentangle <- entanglement(bothdends)
  #cophcors <- c(cophcors,curcor[1,2])
  #entangles <- c(entangles,curentangle)
  tanglegram(bothdends, left_dendo_mar=c(4,4,4,8),right_dendo_mar=c(4,8,4,4),
             main_left = 'DNA', main_right= 'Pathway',
             common_subtrees_color_lines = FALSE)
  mtext(side=3, text=paste0(pat, ': Correlation = ',signif(curcor[1,2],3)),line=0,font=2,cex=0.75)
  dev.off()
  df[df$Patient == pat, 'Cophenetic_corr'] <- signif(curcor[1,2],3)

  } 
}
write.csv(df, file = paste0(outdir, 'cophenetic_correlation_method-', method, '.txt'), quote=F, row.names=F)



#.....................................cytokine genes.......................
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

cytokine_genes <- getBM(
  attributes = c("hgnc_symbol", "description", "chromosome_name", "start_position", "end_position"),
  filters = "go", 
  values = "GO:0005125",  # GO term for cytokine activity
  mart = ensembl
)
head(cytokine_genes)

write.csv(cytokine_genes, "cytokine_genes.csv", row.names = FALSE)


#.....................................COX analysis.......................
#survival analysis
#install.packages(c("survival", "survminer"))
library("survival")
library("survminer")

setwd('./')
data<-read.table('data.txt',header=T,row.names = 1,check.names = F)
#data.txt
#patient_id   PFI    PFI_time    geneA_mutation_status    geneB_mutation_status
#P1            1       100              0                         1
#P2            0        50              1                         1

gene_list<-colnames(data)[5:ncol(data)]

#pdf('PFI-survival-EEF1E1.pdf')
i=4

for (gene in gene_list){
  i=i+1
  cox_model <- coxph(Surv(OS_time, OS) ~ data[,i], data = data)
  summary_cox <- summary(cox_model)
  exp_coef <- summary_cox$coefficients[, "exp(coef)"]
  lower_95 <- summary_cox$conf.int[, "lower .95"]
  upper_95 <- summary_cox$conf.int[, "upper .95"]
  p_value <- summary_cox$coefficients[, "Pr(>|z|)"]
  results <- data.frame(
    SCNA=gene,
    exp_coef = exp_coef,
    lower_95 = lower_95,
    upper_95 = upper_95,
    p_value = p_value
  )
  print(results)
}


#.....................................MCPcounter.......................
#MCPcounter
setwd('./')
library(devtools)
devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")
library(MCPcounter)
expr_matrix <- read.delim('TPM.txt', row.names = 1, header=T,sep = '\t', check.names = FALSE)
results <- MCPcounter.estimate(expr_matrix, featuresType = "HUGO_symbols")
print(results)

#.....................................TIDE.......................
#generate TIDE score
getwd()
setwd('./')
data<-read.table('TPM.txt',header=T,row.names = 1,check.names=FALSE);head(data)

#将表达矩阵均值标准化（行表达量-行均值）
normalize <- t(apply(data, 1, function(x)x-(mean(x))))
write.table(data.frame("gene symbol"=rownames(normalize),normalize),file="TIDE_input-tpm-nor.txt",sep = "\t",quote = F,row.names = F)
#结果可视化
setwd('./')
# 读入结果表：  
result <- read.csv('TIDE-using-TPM-nor.csv')  
colnames(result)

my_comparisons<-list(c('P01_A','P01_BC'))
# 小提琴图展示结果：  
# 1.TIDE小提琴图：  
#TIDE评分越高，肿瘤免疫逃逸的可能性越高，免疫检查点抑制剂疗效越差
p1 <- ggviolin(result, x = 'tumor', y = 'TIDE', fill = 'patient',  
               #palette = c("#2E9FDF", "#E7B800"),  
               add = 'boxplot', add.params = list(fill = "white")) +  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size=0.5, tip.length = 0.02, method = 'wilcox.test')+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
p1  
# 2.Dysfunction小提琴图：
# dysfunction score的计算原理：免疫失调作用的基因拥有更高的权重，再乘以表达量  
p2 <- ggviolin(result, x = 'tumor', y = 'Dysfunction', fill = 'patient',  
               #palette = c("#2E9FDF", "#E7B800"),  
               add = 'boxplot', add.params = list(fill = "white")) +  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size=0.5, tip.length = 0.02, method = 'wilcox.test')+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
p2 
# Exclusion小提琴图：  
# exclusion score是由免疫排斥的基因拥有更高的权重，再乘以表达量得到
p3 <- ggviolin(result, x = 'tumor', y = 'Exclusion', fill = 'patient',  
               #palette = c("#2E9FDF", "#E7B800"),  
               add = 'boxplot', add.params = list(fill = "white")) +  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size=0.5, tip.length = 0.02, method = 'wilcox.test')+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
p3 
p <- p1 + p2 + p3 
p



#.....................................DESEQ2.......................
setwd('./')
library(DESeq2)
library(dplyr)
dat <- read.delim('deseq2-input-only-TNBC.txt', row.names = 1, header=T,sep = '\t', check.names = FALSE)
## QC
### expr_gene_num > 5000 TNBC cel clusters
print(paste('before QC, cell cluster: ',ncol(dat)))
print(paste('before QC, gene: ',length(rownames(dat))))
gene_num <- apply(dat, 2, function(x) { length(x[x>0]) })
gt_5k_sam <- gene_num[gene_num > 5000] %>% names()
dat <- dat[, gt_5k_sam] %>% as.data.frame

### Find which genes are >=5 read count in >=5%  of tumour samples
perc5 <- ceiling(ncol(dat)*.05)
atl10 <- apply(dat,1,function(x) sum(x>=20))
dat$gene <- rownames(dat)
filens <- dat$gene[which(atl10>=perc5)]
dat <- dat[filens,]
dat <- dat[,1:(ncol(dat)-1)]
print(paste('after QC, cell cluster: ',ncol(dat)))
print(paste('after QC, gene: ',length(rownames(dat))))

coldata <- read.delim('deseq2-coldata-only-TNBC.txt',row.names = 1, header=T, sep = '\t', check.names = FALSE)
coldata$condition<-factor(coldata$condition)

dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~condition)

dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'high_cor', 'low_cor'))
res
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'high_low_cor.DESeq2-only-TNBC.txt', col.names = NA, sep = '\t', quote = FALSE)
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res1[which(res1$log2FoldChange > 1.5 & res1$padj < 0.01),'sig'] <- 'up'
res1[which(res1$log2FoldChange < -1.5 & res1$padj < 0.01),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1.5 | res1$padj >= 0.01),'sig'] <- 'none'

res1_select <- subset(res1, sig %in% c('up', 'down'))
res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')
table(res1$sig)
write.table(res1_up, file = 'high_low_cor.DESeq2.up-cutoff15-sig001-only-TNBC.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'high_low_cor.DESeq2.down-cutoff15-sig001-only-TNBC.txt', sep = '\t', col.names = NA, quote = FALSE)

#heatmap of all DEG
library(pheatmap)
plot_df<-dat[rownames(res1_select),]
plot_df<-log2(plot_df+1)
plot_annotation_row<-as.data.frame(res1_select$sig,rownames(res1_select))
plot_annotation_col<-coldata
plot_annotation_col$patient<-sapply(strsplit(rownames(plot_annotation_col),"_"),function(x) x[1])
bk<-seq(-1.5,1.5,by=0.1)
pheatmap(plot_df,
         scale = "row",main="TPM of all DEG in high_corr vs. low_corr",
         cluster_rows = T,show_colnames = F,show_rownames = F,
         cluster_cols = F,border_color = NA,annotation_row = plot_annotation_row,
         annotation_col=plot_annotation_col,
         color = c(colorRampPalette(colors = c("#137AB5","white"))(length(bk)/2),colorRampPalette(colors = c("white","#CC0724"))(length(bk)/2)),
         legend_breaks=seq(-1.5,1.5,0.5),breaks=bk
)

#heatmap of top200 DEG
res1 <- res1[order(res1$log2FoldChange,res1$padj,  decreasing = c(TRUE,FALSE)), ]
res1_up=subset(res1, sig %in% c('up'))[1:200,]
res1_down=subset(res1, sig %in% c('down'))[1:200,]
res1_top200<-rbind(rownames(res1_up),rownames(res1_down))
res1_top200_select<-res1_select[res1_top200,]

write.table(res1_up, file = 'high_low_cor.DESeq2.up-cutoff2-sig001-only-TNBC-top200.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'high_low_cor.DESeq2.down-cutoff2-sig001-only-TNBC-top200.txt', sep = '\t', col.names = NA, quote = FALSE)

plot_df<-dat[res1_top200,]
plot_df<-log2(plot_df+1)
plot_annotation_row<-as.data.frame(res1_top200_select$sig,rownames(res1_top200_select))
plot_annotation_col<-coldata
plot_annotation_col$patient<-sapply(strsplit(rownames(plot_annotation_col),"_"),function(x) x[1])
bk<-seq(-1.5,1.5,by=0.1)
pheatmap(plot_df,
         scale = "row",main="TPM of top200 DEG in high_corr vs. low_corr",
         cluster_rows = T,show_colnames = F,show_rownames = F,
         cluster_cols = F,border_color = NA,annotation_row = plot_annotation_row,
         annotation_col=plot_annotation_col,
         color = c(colorRampPalette(colors = c("#137AB5","white"))(length(bk)/2),colorRampPalette(colors = c("white","#CC0724"))(length(bk)/2)),
         legend_breaks=seq(-1.5,1.5,0.5),breaks=bk
)


library(ggplot2)
library(ggrepel)
genes_to_label<-read.table('./cytokine-list.txt',header=T,stringsAsFactors = FALSE)
genes_to_label <- genes_to_label$Gene
res1$gene<-rownames(res1)
p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  
  scale_color_manual(values = c('#C2395A', 'gray', '#1f78b4'), limits = c('up', 'none', 'down')) +  
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'high_cor vs low_cor', color = '') +  
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 3, color = 'black') +  
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-12, 12) + ylim(0, 35)+  
  geom_text_repel(
  data = subset(res1, gene %in% genes_to_label & res1$sig!='none'),
  aes(label = gene),
  size = 3,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  max.overlaps = 50
)
pdf('high_low_cor.DESeq2.sig-cutoff15-sig001-only-TNBC.pdf')
p
dev.off()

#fgsea
res1<-res1[order(res1$log2FoldChange,decreasing = T),]
id<-res1$log2FoldChange
names(id)<-rownames(res1)

library(fgsea)
Hallmark_cancer_list <- readRDS('./Hallmark_cancer_list.Rds')
names(Hallmark_cancer_list) <- gsub('HALLMARK_(\\S+)','\\1',names(Hallmark_cancer_list))
# Remove spermatogenesis, myogenesis and pancreas_beta_cells
Hallmark_cancer_list$SPERMATOGENESIS <- NULL
Hallmark_cancer_list$MYOGENESIS <- NULL
Hallmark_cancer_list$PANCREAS_BETA_CELLS <- NULL
# Rename COMPLEMENT
Hallmark_cancer_list$COMPLEMENT_INNATE_IMMUNE_SYSTEM <- Hallmark_cancer_list$COMPLEMENT
Hallmark_cancer_list$COMPLEMENT <- NULL

# Include WNT_signalling
geneinfo <- data.table::fread('./complete_gene_info.txt.gz')
Hallmark_cancer_list$WNT_SIGNALING <- read.table('./wnt_signalling_entrez.txt')[,1] %>%intersect(., geneinfo$Entrez) %>%plyr::mapvalues(from = geneinfo$Entrez, to = geneinfo$Symbol)
#write.csv(Hallmark_cancer_list,file='Hallmark_cancer_list_used-this-study.csv')

#unique hallmark gene list
Hallmark_cancer_list<- lapply(Hallmark_cancer_list, unique)

fgseaRes<-fgsea(pathways=Hallmark_cancer_list,
                stats=id)
sig<-fgseaRes[fgseaRes$padj<0.05,]
sig<-sig[order(sig$NES,decreasing=T),]

plotGseaTable(Hallmark_cancer_list[sig$pathway],id,fgseaRes,gseaParam=0.5)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GSEABase")
library(clusterProfiler)
library(GSEABase)
library(stats)
Hallmark_cancer_list<-read.gmt('./Hallmark_cancer_list-unique-gene.gmt')
geneList <- res1$log2FoldChange
names(geneList) <- toupper(rownames(res1))
geneList <- sort(geneList,decreasing = T)

gseaRes<-GSEA(TERM2GENE =Hallmark_cancer_list,
                geneList=geneList,
              pAdjustMethod = 'BH'
              )
write.table(gseaRes, 'gsea_human_hallmark-only-TNBC.txt', sep = '\t', row.names = FALSE, quote = FALSE)


g1<-as.data.frame(gseaRes)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

num<-g1[,c(1,11)]

#enrichplot
library(ggsci)
col_gsea1<-pal_simpsons()(16)

num1=1
gseaplot2(gseaRes,geneSetID = rownames(g1)[1:num1],
          title = "",
          color = col_gsea1[1:num1],
          base_size = 14,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"#line or dot
)

library(enrichplot)
library(ggplot2)


ridgeplot(gseaRes,
          showCategory = 30,
          fill = "pvalue", 
          core_enrichment = TRUE,
          label_format = 30,
          orderBy = "NES",
          decreasing = FALSE
) +theme(axis.text.y = element_text(size=8))