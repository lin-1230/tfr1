######################################
# @date: 2025-04-23
# 这个脚本是为了对Tfr1的单细胞数据进行差异分析和kmeans分析
######################################

#### 在服务器环境中运行
# conda activate Tfr1
# R

## 设置当前工作目录
setwd("/media/ssd/sdb1/data/ljh/TFR1/")


resultDir <- 'result/2025-04-09-zhuShi/'
dir.create(resultDir, recursive = TRUE)
RDSDir <- 'RDS/'
dir.create(RDSDir, recursive = TRUE)


## 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(qs)
library(dplyr)


## 引入配置文件
source("code/config_seurat.R")
## 引入自定义函数
source("code/用到的代码.R")


## 读取数据
SeuratOb <- qread(paste0(RDSDir,'SeuratOb_UMAP_CR_zhuShi.qs'))

table(SeuratOb$cellType)


###### 差异分析 ######


library(parallel)
library(foreach)
library(doParallel)

## 设置并行运行环境
numCores <- detectCores() - 30
registerDoParallel(numCores)

## 将SeuratOb$cellType 中的/替换成_
SeuratOb$cellType <- gsub("/", "_", SeuratOb$cellType)
cellTypes <- unique(SeuratOb$cellType)



DEfunction <- function(SeuratOb,cellClass,resDir,group1,group2,logfc.threshold = 1,p_val.threshold = 0.05,test.use = 'wilcox'){
  ## 参数说明
  ## SeuratOb: Seurat对象
  ## cellClass: 细胞类型,如果是ALL，则使用所有细胞类型
  ## resDir: 结果目录
  ## group1: 分组1，实验组
  ## group2: 分组2，对照组
  ## logfc.threshold: 差异基因的logfc阈值
  ## p_val.threshold: 差异基因的p值阈值
  ## test.use: 差异基因的测试方法，默认为wilcox

  ## 取出细胞类型的细胞
  if (cellClass  == 'ALL'){
   SeuratOb_use <- SeuratOb
  }
  else{
    SeuratOb_use <- subset(SeuratOb, subset = cellType == cellClass)
  }
  
  print(dim(SeuratOb_use))
  ## 重新定义分组
  SeuratOb_use$group2 <- as.character(SeuratOb_use$group)
  ## KO2、KO3、KO4 合并为KO
  SeuratOb_use$group2[SeuratOb_use$group2 %in% c('KO2','KO3','KO4')] <- 'KO'
  ## WT1、WT2合并为WT
  SeuratOb_use$group2[SeuratOb_use$group2 %in% c('WT1','WT2')] <- 'WT'
  
  Idents(SeuratOb_use) <- SeuratOb_use$group2
  ## 运行FindMarkers
  markers <- FindMarkers(SeuratOb_use, ident.1 = group1, ident.2 = group2,min.pct = 0.1,only.pos = F,logfc.threshold = 0,test.use = test.use)
  markers <- markers %>%
    filter(p_val < p_val.threshold)
  
  write.csv(markers,paste0(resDir,cellClass,'_markers.csv'))

  ## 根据阈值筛选显著的基因
  ## 筛选出avg_log2FC > logfc.threshold 且 p_val < p_val.threshold 的上调基因
  markers_use <- markers %>%
    filter(avg_log2FC > logfc.threshold, p_val < p_val.threshold)
  ## 保存
  write.csv(markers_use,paste0(resDir,cellClass,'_markers_up.csv'))

  ## 筛选出avg_log2FC < -logfc.threshold 且 p_val < p_val.threshold 的下调基因
  markers_use <- markers %>%
    filter(avg_log2FC < -logfc.threshold, p_val < p_val.threshold)
  ## 保存
  write.csv(markers_use,paste0(resDir,cellClass,'_markers_down.csv'))
  
  return(markers)
}






## 挑选出每个细胞类型的DEG结果
results <- foreach(cellType = cellTypes, .combine = 'c') %dopar% {
  umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
  cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/',cellType,'/')
  dir.create(cellDir, recursive = TRUE)
  
  DEfunction(SeuratOb = SeuratOb,
             cellClass = cellType,
             resDir = cellDir,
             group1 = 'KO',
             group2 = 'WT',
             logfc.threshold = 1,
             p_val.threshold = 0.05,
             test.use = 'wilcox')
}

## 跑出ALL的DEG结果
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/ALL/')
dir.create(cellDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb,
             cellClass = 'ALL',
             resDir = cellDir,
             group1 = 'KO',
             group2 = 'WT',
             logfc.threshold = 0,
             p_val.threshold = 0.05,
             test.use = 'wilcox')


## 挑选出KO_FAC vs KO ALL的DEG结果
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/ALL_KO_FAC_vs_KO/')
dir.create(cellDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb,
             cellClass = 'ALL',
             resDir = cellDir,
             group1 = 'KO_FAC',
             group2 = 'KO',
             logfc.threshold = 0,
             p_val.threshold = 0.05,
             test.use = 'wilcox')



## 先取出部分基因
k_40_res <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/kmeans_40/kmeans_result.csv')

## 提取出cluster是1 5 8 9 10 20 30 37的基因
k_40_res_use <- k_40_res %>%
  filter(cluster %in% c(1,5,8,9,10,20,30,37))

head(k_40_res_use)
use_genes <- k_40_res_use$rn

## 提取出先上调后下降的基因
SeuratOb_use <- subset(SeuratOb,features = use_genes)
dim(SeuratOb_use)
head(SeuratOb_use)

## 重新定义分组
SeuratOb_use$group2 <- as.character(SeuratOb_use$group)
## KO2、KO3、KO4 合并为KO
SeuratOb_use$group2[SeuratOb_use$group2 %in% c('KO2','KO3','KO4')] <- 'KO'
## WT1、WT2合并为WT
SeuratOb_use$group2[SeuratOb_use$group2 %in% c('WT1','WT2')] <- 'WT'


## 获取按照分组平均后的表达数据
data_use <- AverageExpression(SeuratOb_use, group.by = 'group2', slot = 'data')
data_use <- data_use[[1]]
data_use <- as.data.frame(data_use)
head(data_use)

summary(as.numeric(data_use$KO))
which(data_use$KO > 100)
which(data_use$KO < 100)
summary(as.numeric(data_use$KO))

quantile(data_use$KO,seq(0,1,0.1))




## 跑差异基因（全部细胞类型，基因只用先上调后下降的基因）
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_up/')
dir.create(cellDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb_use,
             cellClass = 'ALL',
             resDir = cellDir,
             group1 = 'KO',
             group2 = 'WT',
             logfc.threshold = 0,
             p_val.threshold = 0.05,
             test.use = 'wilcox')


## 挑选出KO_FAC vs KO ALL的DEG结果
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/ALL_KO_FAC_vs_KO_use_genes_up/')
dir.create(cellDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb_use,
             cellClass = 'ALL',
             resDir = cellDir,
             group1 = 'KO_FAC',
             group2 = 'KO',
             logfc.threshold = 0,
             p_val.threshold = 0.05,
             test.use = 'wilcox')




umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/ALL/')
dir.create(cellDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb,
             cellClass = 'ALL',
             resDir = cellDir,
             group1 = 'KO',
             group2 = 'WT',
             logfc.threshold = 0,
             p_val.threshold = 0.05,
             test.use = 'wilcox')


## 挑选出KO_FAC vs KO ALL的DEG结果
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/ALL_KO_FAC_vs_KO/')
dir.create(cellDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb,
             cellClass = 'ALL',
             resDir = cellDir,
             group1 = 'KO_FAC',
             group2 = 'KO',
             logfc.threshold = 0,
             p_val.threshold = 0.05,
             test.use = 'wilcox')



## 读取差异结果
## 读取KO vs WT的差异结果
DEG_KO_WT <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_up/ALL_markers.csv')
DEG_KO_FAC_KO <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_KO_FAC_vs_KO_use_genes_up/ALL_markers.csv')

## 分别筛选出avg_log2FC > 1的基因
# DEG_KO_WT <- DEG_KO_WT %>%
#   filter(avg_log2FC > 1)
# DEG_KO_FAC_KO <- DEG_KO_FAC_KO %>%
#   filter(avg_log2FC < -1)


## 分别筛选出avg_log2FC > 0.58的基因
# DEG_KO_WT <- DEG_KO_WT %>%
#   filter(avg_log2FC > 0.58)
# DEG_KO_FAC_KO <- DEG_KO_FAC_KO %>%
#   filter(avg_log2FC < -0.58)





## 取交集
DEG <- intersect(DEG_KO_WT$V1,DEG_KO_FAC_KO$V1)

## 保存
write.csv(DEG,'/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_up/DEG_KO_WT_vs_KO_FAC_KO_FC_1.5.csv')


S



## 处理62个基因的富集结果
inputDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/metascape2/DEG_KO_WT_vs_KO_FAC_KO_FC_1.5/'
resDir <- inputDir
plotMetascapeBar(inputDir,resDir,goColor,keggColor)





##############先找到交集的差异基因########################




max(DEG_KO_WT$avg_log2FC)
max(DEG_KO_WT$pct.1)
max(DEG_KO_WT$pct.2)
summary(DEG_KO_WT$pct.1)
summary(DEG_KO_WT$pct.2)



# 停止并行集群
# stopImplicitCluster()



cellType <- 'HSC'
resultDir <-'result/2025-04-09-zhuShi/'
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/',cellType,'/')
 DEfunction(SeuratOb = SeuratOb,
             cellClas = 'HSC',
             resDir = cellDir,
             group1 = 'KO',
             group2 = 'WT',
             logfc.threshold = 1,
             p_val.threshold = 0.05,
             test.use = 'wilcox')

cellType <- 'MPP2_3'
resultDir <-'result/2025-04-09-zhuShi/'
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/',cellType,'/')
 DEfunction(SeuratOb = SeuratOb,
             cellType = 'MPP2_3',
             resDir = cellDir,
             group1 = 'KO',
             group2 = 'WT',
             logfc.threshold = 1,
             p_val.threshold = 0.05,
             test.use = 'wilcox')



## 自定义筛选DEG结果的函数
filter_DEG_results <- function(DEG_results_dir,DEG_results_name, logfc.threshold, p_val.threshold) {
  ## 读取DEG结果
  library(dplyr)

  DEG_results <- read.csv(paste0(DEG_results_dir,DEG_results_name,'.csv'), row.names = 1)
  ## 筛选上调基因
  DEG_results_up <- DEG_results %>%
    filter(avg_log2FC > logfc.threshold, p_val < p_val.threshold)
  ## 筛选下调基因
  DEG_results_down <- DEG_results %>%
    filter(avg_log2FC < -logfc.threshold, p_val < p_val.threshold)

  resDir <- paste0(DEG_results_dir,'筛选结果_',logfc.threshold,'_',p_val.threshold,'/')
  dir.create(resDir, recursive = TRUE)
  ## 保存
  write.csv(DEG_results_up, paste0(resDir,DEG_results_name,'_up.csv'))
  write.csv(DEG_results_down, paste0(resDir,DEG_results_name,'_down.csv'))
}


DEG_results_dir <- '/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/'
## B细胞
B_dir <- paste0(DEG_results_dir,'B/')
DEG_results_name <- 'B_markers'
logfc.threshold <- 0.5
p_val.threshold <- 0.05
filter_DEG_results(B_dir,DEG_results_name, logfc.threshold, p_val.threshold)


## MPP4
MPP4_dir <- paste0(DEG_results_dir,'MPP4/')
DEG_results_name <- 'MPP4_markers'
logfc.threshold <- 0.5
p_val.threshold <- 0.05
filter_DEG_results(MPP4_dir,DEG_results_name, logfc.threshold, p_val.threshold)

## MPP1
MPP1_dir <- paste0(DEG_results_dir,'MPP1/')
DEG_results_name <- 'MPP1_markers'
logfc.threshold <- 0.5
p_val.threshold <- 0.05
filter_DEG_results(MPP1_dir,DEG_results_name, logfc.threshold, p_val.threshold)

## HSC
HSC_dir <- paste0(DEG_results_dir,'HSC/')
DEG_results_name <- 'HSC_markers'
logfc.threshold <- 0.5
p_val.threshold <- 0.05
filter_DEG_results(HSC_dir,DEG_results_name, logfc.threshold, p_val.threshold)

## MPP2_3
MPP2_3_dir <- paste0(DEG_results_dir,'MPP2_3/')
DEG_results_name <- 'MPP2_3_markers'
logfc.threshold <- 0.5
p_val.threshold <- 0.05
filter_DEG_results(MPP2_3_dir,DEG_results_name, logfc.threshold, p_val.threshold)

## cycling_MPP2_3
cycling_MPP2_3_dir <- paste0(DEG_results_dir,'cycling_MPP2_3/')
DEG_results_name <- 'cycling_MPP2_3_markers'
logfc.threshold <- 0.5
p_val.threshold <- 0.05
filter_DEG_results(cycling_MPP2_3_dir,DEG_results_name, logfc.threshold, p_val.threshold)


################### 绘制整理的metascape通路富集结果

## 本地环境
## conda activate seurat
## R


plotMetascape <- function(upDir,downDir,resDir,color) {
  
  ## 读取metascape结果
  library(readxl)
  library(data.table)
  metascape_results_up <- readxl::read_xlsx(paste0(upDir,'metascape_result.xlsx'),sheet='Enrichment')
  metascape_results_down <- readxl::read_xlsx(paste0(downDir,'metascape_result.xlsx'),sheet='Enrichment')

  metascape_results_up <- as.data.frame(metascape_results_up)
  metascape_results_down <- as.data.frame(metascape_results_down)

  ## 筛选出GroupID带有Summary的行
  metascape_results_up <- metascape_results_up[grep('Summary', metascape_results_up$GroupID),]
  metascape_results_down <- metascape_results_down[grep('Summary', metascape_results_down$GroupID),]

  metascape_results_up$color <- 'UP'
  metascape_results_down$color <- 'DOWN'
  ## 排序
  metascape_results_up <- metascape_results_up[order(metascape_results_up$`LogP`),]
  metascape_results_down <- metascape_results_down[order(metascape_results_down$`LogP`),]
  
  ## 合并两个dataframe
  metascape_results <- rbind(metascape_results_up,metascape_results_down)

  metascape_results$Description <- factor(metascape_results$Description, levels = metascape_results$Description)

  ## 绘制bart图，横坐标是-log10(p-value)，纵坐标是Description，颜色是color
  library(ggplot2)

  ## 绘制
  pdf(resDir, width = pdf_width*2, height = pdf_height)
  p <- ggplot(metascape_results, aes(x = -LogP, y = Description, fill = color)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw() +
    labs(x = '-log10(p-value)', y = 'Description') +
    ## 定义颜色，UP是红色，DOWN是蓝色
    scale_fill_manual(values = color,
                      name = 'Change',
                      labels = c('UP', 'DOWN'))  & themeSet
  print(p)
  dev.off() 

}



## HSC
upDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/HSC_up/'
downDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/HSC_down/'
resDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/HSC_metascape.pdf'
color <- c('red','blue')
plotMetascape(upDir,downDir,resDir,color)

## MPP1
upDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP1_up/'
downDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP1_down/' 
resDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP1_metascape.pdf'
color <- c('red','blue')
plotMetascape(upDir,downDir,resDir,color)
## MPP2_3
upDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP2_3_up/'
downDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP2_3_down/' 
resDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP2_3_metascape.pdf'
color <- c('red','blue')
plotMetascape(upDir,downDir,resDir,color)

## cycling_MPP2_3
upDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/cycling_MPP2_3_up/'
downDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/cycling_MPP2_3_down/'
resDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/cycling_MPP2_3_metascape.pdf'
color <- c('red','blue')
plotMetascape(upDir,downDir,resDir,color)

## MPP4
upDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP4_up/'
downDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP4_down/'
resDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/2025-04-21-Tfr1表达/DEG/metascape/MPP4_metascape.pdf'
color <- c('red','blue')
plotMetascape(upDir,downDir,resDir,color)


### kmeans分析############

## 绘制热图
library(pheatmap)
library(data.table)


## 新建一个group2列
## WT1和WT2的合并成WT
## KO2、KO3和KO4的合并成KO
SeuratOb$group2 <- as.character(SeuratOb$group)
SeuratOb$group2[SeuratOb$group2 %in% c('WT1','WT2')] <- 'WT'
SeuratOb$group2[SeuratOb$group2 %in% c('KO2','KO3','KO4')] <- 'KO'
table(SeuratOb$group2)

## 获取按照分组平均后的表达数据
data <- AverageExpression(SeuratOb, group.by = 'group2', slot = 'data')

data <- data[[1]]
## 调整data列的顺序：WT，KO，KO_FAC，WT_FAC
data <- data[,c('WT','KO','KO_FAC','WT_FAC')]

head(data)



## 寻找 KO > WT 且 KOFAC < KO 的基因
diff_genes_010 <- data.frame(
  gene = rownames(data),
  KO_vs_WT = data[,"KO"] - data[,"WT"],
  KOFAC_vs_KO = data[,"KO_FAC"] - data[,"KO"]
) %>% 
  filter(KO_vs_WT > 0, KOFAC_vs_KO < 0)  # 调整阈值

## 绘制热图（只绘制有差异的基因）
pdf(paste0(resultDir,'UMAP-harmony_CR/热图_小大小.pdf'), width = pdf_width, height = pdf_height)
pheatmap(data[diff_genes_010$gene,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 10,
         fontsize_col = 10,
         #legend_breaks = c(0, 0.5, 1),
         #legend_labels = c("Low", "Medium", "High"),
         main = ""
)
dev.off()



## 绘制热图（全部基因）
pdf(paste0(resultDir,'UMAP-harmony_CR/热图_all.pdf'), width = pdf_width, height = pdf_height)
pheatmap(data,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 10,
         fontsize_col = 10,
         #legend_breaks = c(0, 0.5, 1),
         #legend_labels = c("Low", "Medium", "High"),
         main = ""
)
dev.off()

## 绘制折线图
pdf(paste0(resultDir,'UMAP-harmony_CR/折线图.pdf'), width = pdf_width, height = pdf_height)
plot(diff_genes,
     type = "l",
     lty = 1,
     col = "black",,
     type = "l",
     lty = 1,
     col = "black",
     xlab = "Samples",
     ylab = "Expression",
     main = "Gene Expression Over Samples"
) 
dev.off()


library(ggplot2)
pdf(paste0(resultDir,'UMAP-harmony_CR/折线图.pdf'), width = pdf_width, height = pdf_height)
ggplot(data, aes(x = , y = KO_vs_WT)) +
  geom_line() +
  labs(x = "Gene", y = "Expression", title = "Gene Expression Changes")
dev.off()



## 进行kmeans聚类，然后取出其中每一类的基因
set.seed(666)
kmeans_result <- kmeans(data, centers = 12, nstart = 25, iter.max = 100)
kmeans_result$cluster

## 查看每个聚类的基因数量
table(kmeans_result$cluster)
head(data)
kmeans_result$cluster[1:5]

## 绘制每个聚类的基因表达的折线图，使用ggplot2
## 首先，将聚类结果添加到data数据框中
data <-  as.data.table(data,keep.rownames = TRUE)
data$cluster <- kmeans_result$cluster
head(data)



## 绘制分组折线图，每个线是一个基因，横坐标是WT，KO，KO_FAC，WT_FAC，纵坐标是表达值

## 首先改变data的格式适合ggplot2
library(tidyr)
data_gather <- data %>%
  gather(key = "group2", value = "value", -rn, -cluster)
table(data$group2)
table(data$cluster)

data$group2 <- factor(data$group2, levels = c("WT", "KO", "KO_FAC", "WT_FAC"))


for (cluster in unique(data$cluster)) {
  ## 筛选出cluster为cluster的行
  data_cluster <- data[data$cluster == cluster,]
  ## 绘制分组折线图，每个线是一个基因，横坐标是WT，KO，KO_FAC，WT_FAC，纵坐标是表达值
  pdf(paste0(resultDir,'UMAP-harmony_CR/分组折线图_cluster_',cluster,'.pdf'), width = pdf_width, height = pdf_height)
  p <- ggplot(data_cluster, aes(x = group2, y = value, group = rn)) +
    geom_line() +
    labs(x = "Group", y = "Expression", title = "Gene Expression Changes") & themeSet
  print(p) 
  dev.off()
}







## 获取按照分组平均后的表达数据
data <- AverageExpression(SeuratOb, group.by = 'group2', slot = 'data')

data <- data[[1]]
## 调整data列的顺序：WT，KO，KO_FAC，WT_FAC
data <- data[,c('WT','KO','KO_FAC','WT_FAC')]

## 对data进行行scale
data <- t(scale(t(data)))
head(data)



plotKmeansLine <- function(data,resultDir,k=20,seed=666,nstart=25,iter.max=100) {
  library(ggplot2)
  library(tidyr)
  library(data.table)

  ## 进行kmeans聚类，然后取出其中每一类的基因
  set.seed(seed)
  kmeans_result <- kmeans(data, centers = k, nstart = nstart, iter.max = iter.max)

  ## 查看每个聚类的基因数量
  print(table(kmeans_result$cluster))
  data <-  as.data.table(data,keep.rownames = TRUE)
  data$cluster <- kmeans_result$cluster
  data_gather <- data %>%
    gather(key = "group2", value = "value", -rn, -cluster)

  data_gather$group2 <- factor(data_gather$group2, levels = c("WT", "KO", "KO_FAC", "WT_FAC"))
  for (cluster in unique(data_gather$cluster)) {
  ## 筛选出cluster为cluster的行
  data_cluster <- data_gather[data_gather$cluster == cluster,]
  ## 绘制分组折线图，每个线是一个基因，横坐标是WT，KO，KO_FAC，WT_FAC，纵坐标是表达值
  pdf(paste0(resultDir,'分组折线图_cluster_',cluster,'.pdf'), width = pdf_width, height = pdf_height)
  p <- ggplot(data_cluster, aes(x = group2, y = value, group = rn)) +
    geom_line(color = "steelblue",alpha=0.6) +
    labs(x = "Group", y = "Expression", title = "Gene Expression Changes") & themeSet
  print(p) 
  dev.off()
  }
  write.csv(data,paste0(resultDir,'kmeans_result.csv'))
}

## k=20
k <- 20
resultDir <- paste0('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/','UMAP-harmony_CR/kmeans_',k,'/')
dir.create(resultDir, recursive = TRUE)
plotKmeansLine(data,resultDir,k=k,seed=666,nstart=25,iter.max=100)

## k=10
k <- 10
resultDir <- paste0('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/','UMAP-harmony_CR/kmeans_',k,'/')
dir.create(resultDir, recursive = TRUE)
plotKmeansLine(data,resultDir,k=k,seed=666,nstart=25,iter.max=100)

## k=30 
k <- 30
resultDir <- paste0('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/','UMAP-harmony_CR/kmeans_',k,'/')
dir.create(resultDir, recursive = TRUE)
plotKmeansLine(data,resultDir,k=k,seed=666,nstart=25,iter.max=100)

## k=40
k <- 40
resultDir <- paste0('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/','UMAP-harmony_CR/kmeans_',k,'/')
dir.create(resultDir, recursive = TRUE)
plotKmeansLine(data,resultDir,k=k,seed=666,nstart=25,iter.max=100)






 ## 自定义绘制KEGG的通路富集bar图的函数
plotMetascapeBar <- function(inputDir,resDir,goColor,keggColor){
  
  ## 读取metascape结果
  library(readxl)
  library(data.table)
  metascape_results <- readxl::read_xlsx(paste0(inputDir,'metascape_result.xlsx'),sheet='Enrichment')
  
  metascape_results <- as.data.frame(metascape_results)
  
  ## 去掉出GroupID带有Summary的行
  metascape_results <- metascape_results[!grepl('Summary', metascape_results$GroupID),]


  ## 绘制Go部分
  ## 筛选出Category带有GO的行
  metascape_results_go <- metascape_results[grep('GO', metascape_results$Category),]
  metascape_results_go$color <- metascape_results_go$Category
  
  ## 排序
  metascape_results_go <- metascape_results_go[order(-metascape_results_go$`LogP`),]
  metascape_results_go$Description <- factor(metascape_results_go$Description, levels = metascape_results_go$Description)

  ## 如果通路数量大于80,只保留前80个通路
  if (nrow(metascape_results_go) > 80) {
    metascape_results_go <- metascape_results_go[1:80,]
  }
  

  ## 绘制GO部分的bar图
  ## 上调的bar图
  pdf(paste0(resDir,'GO.pdf'), width = pdf_width*3, height = pdf_height*3)
  p <- ggplot(metascape_results_go, aes(x = -LogP, y = Description, fill = color)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw() +
    labs(x = '-log10(p-value)', y = 'Description') +
    ## 定义颜色，
    scale_fill_manual(values = goColor,
                      name = 'Category',
                      labels = c('Biological Processes', 'Cellular Components','Molecular Functions'))  & themeSet +
                      ## 修改字体大小
                      theme(axis.text = element_text(size = 30),
                            axis.title = element_text(size = 30),
                            legend.text = element_text(size = 30),
                            legend.title = element_text(size = 30))                      
  print(p)
  dev.off()

  
  ## 绘制KEGG部分
  ## 筛选出Category带有KEGG的行
  metascape_results_kegg <- metascape_results[grep('KEGG', metascape_results$Category),]
  ## 排序
  metascape_results_kegg <- metascape_results_kegg[order(-metascape_results_kegg$`LogP`),]
  ## 
  metascape_results_kegg$Description <- factor(metascape_results_kegg$Description, levels = metascape_results_kegg$Description)

  ## 绘制KEGG部分的bar图
  pdf(paste0(resDir,'KEGG.pdf'), width = pdf_width * 1.2, height = pdf_height)
  p <- ggplot(metascape_results_kegg, aes(x = -LogP, y = Description, fill = keggColor)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw() +
    labs(x = '-log10(p-value)', y = 'Description') & themeSet
  print(p)
  dev.off()

}


goColor <- c('#8BB8E8','#FFC470','#7BC8A4')
keggColor <- 'blue'
metaDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/kemans/metascape_40/'
## up 1 5 8 9 10 20 30 37
## down 3 7 11 16 19 22 33 36 39
for (cluster in paste0('C',c(1,5,8,9,10,20,30,37,3,7,11,16,19,22,33,36,39))) {
  inputDir <- paste0(metaDir,cluster,'/')
  resDir <- paste0(metaDir,cluster,'/')
  plotMetascapeBar(inputDir,resDir,goColor,keggColor)
}



goColor <- c('#8BB8E8','#FFC470','#7BC8A4')
metaDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/metascape/'

cellClass <- c('MPP1','MPP2_3','HSC','MPP4','cycling_MPP2_3')
for (cell in cellClass) {
  inputDir <- paste0(metaDir,cell,'_markers_up/')
  resDir <- paste0(metaDir,cell,'_markers_up/')
  plotMetascapeBar(inputDir,resDir,goColor,keggColor)
  # print(inputDir)

  inputDir <- paste0(metaDir,cell,'_markers_down/')
  resDir <- paste0(metaDir,cell,'_markers_down/')
  plotMetascapeBar(inputDir,resDir,goColor,keggColor)

}



## C37 + C8
inputDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/kemans/metascape_40/C8_C37/'
resDir <- inputDir
plotMetascapeBar(inputDir,resDir,goColor,keggColor)

## 7 11 36 
inputDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/kemans/metascape_40/C7_C11_C36/'
resDir <- inputDir
plotMetascapeBar(inputDir,resDir,goColor,keggColor)





## 