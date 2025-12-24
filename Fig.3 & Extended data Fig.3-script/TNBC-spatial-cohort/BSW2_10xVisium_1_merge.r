#.libPaths("/hsfscqjf2/ztron/P22Z10200N0623/P22Z10200N0623/changeDir/lianzhiwei1/software/miniconda/Or/envs/seuratV5/lib/R/library")
library(ggplot2)
library(Seurat);packageVersion("Seurat")
library(dplyr)
library(tidyr)
library(tibble)
#library(ggpubr)
#library(ggsci)
library(pheatmap)
library(harmony)

# 涉及多样本整合，需要提前设置大一点的内存
library(future)
options(future.globals.maxSize = 500000 * 1024^2) # 设置为1GB
plan("multicore", workers = 6) # 根据您的服务器配置调整



data_dir_list <- c("/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sAA1/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sAA2/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sAA6/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sAA8/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sAA9/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sEA1/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sEA2/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sEA5/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sEA6/outs/",
"/hsfscqjf1/ST_CQ/P24Z32300N0033/GT_data/PMID-40675986-ST/10x.visium/sEA7/outs/")
length(data_dir_list)
names(data_dir_list) <- c("sAA1", "sAA2", "sAA6", "sAA8", "sAA9", "sEA1", "sEA2", "sEA5", "sEA6", "sEA7")



# 添加样本对应的分组
samg <- read.table(file = "BSW2_sample_group.txt", header = T, sep = "\t")

all_visium <- list()

for(samid in names(data_dir_list)){
    #方案三：包含图像数据的完整读取
    # 完整读取Visium数据（包含图像）
    visium_data <- Load10X_Spatial(
      data.dir = data_dir_list[[samid]],
      filename = "filtered_feature_bc_matrix.h5",
      assay = "Spatial",
      slice = samid,
      filter.matrix = TRUE,
      image = NULL  # 自动读取图像
    )
    
    # 或者手动添加图像
    image <- Read10X_Image(
      image.dir = file.path(data_dir_list[[samid]], "spatial"),
      image.name = "tissue_lowres_image.png"
    )
    visium_data[[samid]] <- image
    visium_data$sample <- samid
    visium_data$orig.ident <- samid
    visium_data[["race"]] <- samg %>% filter(sampleid == samid) %>% pull(group) %>% unique()

    # 质量控制流程
    # 计算线粒体基因比例
    visium_data[["percent.mt"]] <- PercentageFeatureSet(visium_data, pattern = "^MT-")
    # 计算核糖体基因比例
    visium_data[["percent.ribo"]] <- PercentageFeatureSet(visium_data, pattern = "^RP[SL]")
    # 查看质量指标
    visium_data <- subset(visium_data, subset = nFeature_Spatial > 200 & nCount_Spatial > 500 & percent.mt < 25 & percent.ribo < 25)
    
    # 空间分布可视化
    p1 <- SpatialFeaturePlot(visium_data, features = "nCount_Spatial") + theme(legend.position = "right")
    p2 <- SpatialFeaturePlot(visium_data, features = "nFeature_Spatial") + theme(legend.position = "right")
    ggsave(plot = p1 + p2, filename = paste0("BSW2_SpatialFeature_", samid, ".pdf"))
    
    print(samid)
    
    all_visium[[samid]] <- visium_data
}

object.size(all_visium)



# 使用标准的SCTransform流程
all_visium <- lapply(all_visium, function(x) {
  SCTransform(x, assay = "Spatial", verbose = FALSE)
})

saveRDS(all_visium, file = "BSW2_split_SCTransform.rds")


# 2. 选择整合特征并合并
#features <- SelectIntegrationFeatures(object.list = all_visium, nfeatures = 3000)
# merge samples into one object
# 合并所有样本
merged_seurat <- merge(
  all_visium[[1]],
  y = all_visium[-1],
  add.cell.ids = names(all_visium),
  project = "BSW2_10Xvisium"
)


# 3. 运行PCA
merged_seurat <- SCTransform(merged_seurat, assay = "Spatial", verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, npcs = 50, verbose = FALSE) #features = features, 

####################################### 4. Harmony批次校正
merged_seurat <- RunHarmony(merged_seurat, group.by.vars = "sample", assay.use = "SCT") #max.iter.harmony = 20 可调整迭代次数


# 5. 下游分析
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = c(0.2, 0.4, 0.6, 0.8))

saveRDS(merged_seurat, file = "BSW2_integrate_harmony.rds")

# 检查批次效应校正效果
#p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "sample") 
p2 <- DimPlot(merged_seurat, reduction = "umap", label = TRUE)
pdf("BSW2_integrate_FindClusters_Dimplot.pdf", width = 6, height = 6)
#print(p1)
print(p2)
dev.off()