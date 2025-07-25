######################################
# @date: 2025-04-25
# 这个脚本是为了对Tfr1的单细胞数据进行kmeans分析
# 只对HSC进行分析
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


## 提取出HSC的细胞
SeuratOb <- subset(SeuratOb, subset = cellType == 'HSC')
dim(SeuratOb)
SeuratOb$group2 <- as.character(SeuratOb$group)
## KO2、KO3、KO4 合并为KO
SeuratOb$group2[SeuratOb$group2 %in% c('KO2','KO3','KO4')] <- 'KO'
## WT1、WT2合并为WT
SeuratOb$group2[SeuratOb$group2 %in% c('WT1','WT2')] <- 'WT'




## 获取按照分组平均后的表达数据
data <- AverageExpression(SeuratOb, group.by = 'group2', slot = 'data')

data <- data[[1]]
data <- data[rowSums(data) > 0, ] ## 去除表达值为0的基因

## 调整data列的顺序：WT，KO，KO_FAC，WT_FAC
data <- data[,c('WT','KO','KO_FAC','WT_FAC')]

## 对data进行行scale
data <- t(scale(t(data)))
which(is.nan(data))
head(data)

rowSums(data)

## 自定义函数，用于绘制kmeans聚类的折线图
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

## k=10
k <- 10
resultDir <- paste0('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/','UMAP-harmony_CR/kmeans_HSC_',k,'/')
dir.create(resultDir, recursive = TRUE)
plotKmeansLine(data,resultDir,k=k,seed=666,nstart=25,iter.max=100)

## k=20
k <- 20
resultDir <- paste0('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/','UMAP-harmony_CR/kmeans_HSC_',k,'/')
dir.create(resultDir, recursive = TRUE)
plotKmeansLine(data,resultDir,k=k,seed=666,nstart=25,iter.max=100)

## k=30
k <- 30
resultDir <- paste0('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/','UMAP-harmony_CR/kmeans_HSC_',k,'/')
dir.create(resultDir, recursive = TRUE)
plotKmeansLine(data,resultDir,k=k,seed=666,nstart=25,iter.max=100)

## k=40
k <- 40
resultDir <- paste0('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/','UMAP-harmony_CR/kmeans_HSC_',k,'/')
dir.create(resultDir, recursive = TRUE)
plotKmeansLine(data,resultDir,k=k,seed=666,nstart=25,iter.max=100)


## up 13 15 16 17 21 22 24 34 
## down 6 14 20 26 31 32 33 35 37 38 



## 寻找先下调后上调的基因
k_40_res <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/kmeans_HSC_40/kmeans_result.csv')

## 读取Seurat对象
dim(SeuratOb)


## 提取出cluster是 6 14 20 26 31 32 33 35 37 38 
k_40_res_use <- k_40_res %>%
  filter(cluster %in% c(6,14,20,26,31,32,33,35,37,38))

head(k_40_res_use)
use_genes <- k_40_res_use$rn
dim(k_40_res_use)

## 提取出先上调后下降的基因
SeuratOb_use <- subset(SeuratOb,features = use_genes)
dim(SeuratOb_use)
head(SeuratOb_use)

table(SeuratOb_use$group2 )

# ## 重新定义分组
# SeuratOb_use$group2 <- as.character(SeuratOb_use$group)
# ## KO2、KO3、KO4 合并为KO
# SeuratOb_use$group2[SeuratOb_use$group2 %in% c('KO2','KO3','KO4')] <- 'KO'
# ## WT1、WT2合并为WT
# SeuratOb_use$group2[SeuratOb_use$group2 %in% c('WT1','WT2')] <- 'WT'


## 获取按照分组平均后的表达数据
data_use <- AverageExpression(SeuratOb_use, group.by = 'group2', slot = 'data')
data_use <- data_use[[1]]
data_use <- as.data.frame(data_use)
head(data_use)
dim(data_use)


dim(SeuratOb_use)



## 跑差异基因（全部细胞类型，基因先下降后上调的基因）
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_use_genes_down/')
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
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_KO_FAC_vs_KO_use_genes_down/')
dir.create(cellDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb_use,
             cellClass = 'ALL',
             resDir = cellDir,
             group1 = 'KO_FAC',
             group2 = 'KO',
             logfc.threshold = 0,
             p_val.threshold = 0.05,
             test.use = 'wilcox')



## 读取差异结果
## 读取KO vs WT的差异结果
DEG_KO_WT <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_use_genes_down/ALL_markers.csv')
DEG_KO_FAC_KO <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_KO_FAC_vs_KO_use_genes_down/ALL_markers.csv')

## 筛选出p值小于0.05的基因
DEG_KO_WT <- DEG_KO_WT %>%
  filter(p_val < 0.05)
DEG_KO_FAC_KO <- DEG_KO_FAC_KO %>%
  filter(p_val < 0.05)


max(DEG_KO_WT$avg_log2FC)
min(DEG_KO_WT$avg_log2FC)
max(DEG_KO_FAC_KO$avg_log2FC)
min(DEG_KO_FAC_KO$avg_log2FC)


DEG <- intersect(DEG_KO_WT$V1,DEG_KO_FAC_KO$V1)
length(DEG)

## 提取含有差异基因的结果
data_DEG <- data_use[DEG,]
dim(data_DEG)

head(DEG_KO_WT)


# 分别筛选出avg_log2FC > 1的基因并合并
DEG_KO_WT_1 <- DEG_KO_WT %>%
  filter(avg_log2FC < -1)
DEG_KO_FAC_KO_1 <- DEG_KO_FAC_KO %>%
  filter(avg_log2FC > 1)

DEG_1 <- intersect(DEG_KO_WT_1$V1,DEG_KO_FAC_KO_1$V1)

loc <- which(row.names(data_DEG) %in% DEG_1)
row.names(data_DEG)[loc]


# 分别筛选出avg_log2FC > 0.58的基因并合并
DEG_KO_WT_0p58 <- DEG_KO_WT %>%
  filter(avg_log2FC < -0.58)
DEG_KO_FAC_KO_0p58 <- DEG_KO_FAC_KO %>%
  filter(avg_log2FC > 0.58)


DEG_0p58 <- intersect(DEG_KO_WT_0p58$V1,DEG_KO_FAC_KO_0p58$V1)
loc <- which(row.names(data_DEG) %in% DEG_0p58)
row.names(data_DEG)[loc]

## 保存结果
write.csv(row.names(data_DEG)[loc],'/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_down/DEG_KO_WT_vs_KO_FAC_KO_FC_1.5.csv')



# 分别筛选出avg_log2FC > 0.3的基因并合并
DEG_KO_WT_0p3 <- DEG_KO_WT %>%
  filter(avg_log2FC < -0.3)
DEG_KO_FAC_KO_0p3 <- DEG_KO_FAC_KO %>%
  filter(avg_log2FC >0.3)
DEG_0p3 <- intersect(DEG_KO_WT_0p3$V1,DEG_KO_FAC_KO_0p3$V1)
loc <- which(row.names(data_DEG) %in% DEG_0p3)
row.names(data_DEG)[loc]
write.csv(row.names(data_DEG)[loc],'/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_use_genes_down/DEG_KO_WT_vs_KO_FAC_KO_FC_1.2.csv')




## 寻找先上调后下降的基因
k_40_res <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/kmeans_HSC_40/kmeans_result.csv')

## 读取Seurat对象
dim(SeuratOb)
## 提取出cluster是 13 15 16 17 21 22 24 34
k_40_res_use <- k_40_res %>%
  filter(cluster %in% c(13,15,16,17,21,22,24,34))

head(k_40_res_use)
use_genes <- k_40_res_use$rn
dim(k_40_res_use)
## 提取出先上调后下降的基因
SeuratOb_use <- subset(SeuratOb,features = use_genes)
dim(SeuratOb_use)

## 获取按照分组平均后的表达数据
data_use <- AverageExpression(SeuratOb_use, group.by = 'group2', slot = 'data')
data_use <- data_use[[1]]
data_use <- as.data.frame(data_use)
head(data_use)
dim(data_use)

dim(SeuratOb_use)
## 跑差异基因（全部细胞类型，基因先下降后上调的基因）
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_use_genes_up/')
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
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_KO_FAC_vs_KO_use_genes_up/')
dir.create(cellDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb_use,
             cellClass = 'ALL',
             resDir = cellDir,
             group1 = 'KO_FAC',
             group2 = 'KO',
             logfc.threshold = 0,
             p_val.threshold = 0.05,
             test.use = 'wilcox')


## 读取差异结果
## 读取KO vs WT的差异结果
DEG_KO_WT <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_use_genes_up/ALL_markers.csv')
DEG_KO_FAC_KO <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_KO_FAC_vs_KO_use_genes_up/ALL_markers.csv')
## 筛选出p值小于0.05的基因
DEG_KO_WT <- DEG_KO_WT %>%
  filter(p_val < 0.05)
DEG_KO_FAC_KO <- DEG_KO_FAC_KO %>%
  filter(p_val < 0.05)

max(DEG_KO_WT$avg_log2FC)
min(DEG_KO_WT$avg_log2FC)
max(DEG_KO_FAC_KO$avg_log2FC) 
min(DEG_KO_FAC_KO$avg_log2FC)

DEG <- intersect(DEG_KO_WT$V1,DEG_KO_FAC_KO$V1)
length(DEG)

## 提取含有差异基因的结果
data_DEG <- data_use[DEG,]
dim(data_DEG)
head(DEG_KO_WT)
# 分别筛选出avg_log2FC > 1的基因并合并
DEG_KO_WT_1 <- DEG_KO_WT %>%
  filter(avg_log2FC > 1)
DEG_KO_FAC_KO_1 <- DEG_KO_FAC_KO %>%
  filter(avg_log2FC < -1)
DEG_1 <- intersect(DEG_KO_WT_1$V1,DEG_KO_FAC_KO_1$V1)
loc <- which(row.names(data_DEG) %in% DEG_1)
row.names(data_DEG)[loc]

# 分别筛选出avg_log2FC > 0.58的基因并合并
DEG_KO_WT_0p58 <- DEG_KO_WT %>%
  filter(avg_log2FC > 0.58) 
DEG_KO_FAC_KO_0p58 <- DEG_KO_FAC_KO %>%
  filter(avg_log2FC < -0.58)
DEG_0p58 <- intersect(DEG_KO_WT_0p58$V1,DEG_KO_FAC_KO_0p58$V1)
loc <- which(row.names(data_DEG) %in% DEG_0p58)
row.names(data_DEG)[loc]
# 分别筛选出avg_log2FC > 0.3的基因并合并
DEG_KO_WT_0p3 <- DEG_KO_WT %>%
  filter(avg_log2FC > 0.3)
DEG_KO_FAC_KO_0p3 <- DEG_KO_FAC_KO %>%
  filter(avg_log2FC < -0.3)
DEG_0p3 <- intersect(DEG_KO_WT_0p3$V1,DEG_KO_FAC_KO_0p3$V1)
loc <- which(row.names(data_DEG) %in% DEG_0p3)
row.names(data_DEG)[loc]
## 保存结果
write.csv(row.names(data_DEG)[loc],'/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG_HSC/ALL_use_genes_up/DEG_KO_WT_vs_KO_FAC_KO_FC_1.2.csv')



inputDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG_HSC/metascape/HSC_down/'
resDir <- inputDir
plotMetascapeBar(inputDir,resDir,goColor,keggColor)


inputDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG_HSC/metascape/HSC_up/'
resDir <- inputDir
plotMetascapeBar(inputDir,resDir,goColor,keggColor)


