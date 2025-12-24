library(Seurat)
library(spacexr)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(Matrix.utils)
library(future)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(argparse)
library(ggsci)
library(RColorBrewer)
library(scales)
library(ComplexHeatmap)

# 设置参数
sample = "sAA9"
out = "/ST-TNBC/"
rds = paste0(out, sample, ".rds")
wk = paste0(out, sample, "/")  # 工作目录

# 创建工作目录并设置
if(!dir.exists(wk)){
  dir.create(wk)
  print("工作目录已创建")    
} else{
  print("工作目录已存在")    
}
setwd(wk)

st <- readRDS("/ST-TNBC/sAA9.rds")    
DefaultAssay(st) <- "Spatial" 

# 创建RCTD对象
spatial_counts <- GetAssayData(st, slot = "counts")

coords<-data.frame(x = st@meta.data$x_lowres, y = st@meta.data$y_lowres)
rownames(coords) = rownames(st@meta.data)
names(coords) <- c('x', 'y')

spatial_nUMI <- colSums(spatial_counts)


puck <- SpatialRNA(coords = coords,
                   counts = spatial_counts,
                   nUMI = spatial_nUMI)

counts <- readRDS("/ST-TNBC/ref_counts.rds")
cluster <- read.csv("/ST-TNBC/ref_metadata.csv",row.names=1,header=T)
cluster <- as.factor(cluster$cellType)
names(cluster) <- colnames(counts) 
nUMI <- colSums(counts)

reference <- Reference(counts = counts,
                       cell_types = cluster,
                       nUMI = nUMI)

myRCTD <- create.RCTD(puck, reference, counts_MIN =0,UMI_min=0,max_cores = 20)

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #运行RCTD

save(myRCTD,file= paste0(sample,"_cellbin_myRCTD.Rdata"))


load(file= "/ST-TNBC/sAA9/sAA9_cellbin_myRCTD.Rdata")

results <- myRCTD@results

norm_weights = normalize_weights(results$weights)

saveRDS(norm_weights, file= "cellbin.norm_weights.rds")

cell_type_names <- myRCTD@cell_type_info$info[[2]]

spatialRNA <- myRCTD@spatialRNA
dim(results$results_df)
table(results$results_df$spot_class)
table(results$results_df$first_type)
table(results$results_df$second_type)
head(norm_weights)

spatialRNA <- myRCTD@spatialRNA
cell <- readRDS(file="cellbin.norm_weights.rds")
max_cell <- apply(cell, 1, function(i){
  colnames(cell)[which.max(i)]
})
table(max_cell)

data<- data.frame(
  celltype = max_cell,
  samples = names(max_cell),
  xcoord = spatialRNA@coords[, 1],  
  ycoord = spatialRNA@coords[, 2]
)

fist_cell <- function(x){
  a = order(as.numeric(x),decreasing = T)
  a1 = as.numeric(a[1])
  a2 = as.numeric(a[2])
  
  c = data.frame(
    fist_cell = names(x)[a1],
    fist_weights = x[a1],
    second_cell = names(x)[a2],
    second_weights = x[a2]
  )
  return(c)
}

weights = readRDS("cellbin.norm_weights.rds")
weights = as.data.frame(weights)
fist_weights <- apply(weights,1,fist_cell)
cell_weights = Reduce(rbind,fist_weights)
rownames(cell_weights) = names(fist_weights)
cell_weights$ID = names(fist_weights)

spatialRNA <- myRCTD@spatialRNA
coords <- spatialRNA@coords
cell_weights$x_coord <- coords[cell_weights$ID, 1]
cell_weights$y_coord <- coords[cell_weights$ID, 2]

write.csv(cell_weights,"cell_weights.csv")