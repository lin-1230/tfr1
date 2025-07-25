##############先找到交集的差异基因########################


## 读取差异结果
## 读取KO vs WT的差异结果
DEG_KO_WT <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_up/ALL_markers.csv')
DEG_KO_FAC_KO <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_KO_FAC_vs_KO_use_genes_up/ALL_markers.csv')

## 筛选出p值小于0.05的基因
DEG_KO_WT <- DEG_KO_WT[DEG_KO_WT$p_val_adj < 0.05,]
DEG_KO_FAC_KO <- DEG_KO_FAC_KO[DEG_KO_FAC_KO$p_val_adj < 0.05,]



## 取交集
DEG <- intersect(DEG_KO_WT$V1,DEG_KO_FAC_KO$V1)

## 保存
write.csv(DEG,'/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_up/DEG_KO_WT_vs_KO_FAC_KO_FC_1.5.csv')



## 提取含有差异基因的结果
data_DEG <- data_use[DEG,] 
dim(data_DEG)

quantile(data_DEG$KO,seq(0,1,0.1))

data_DEG <- data_DEG[data_DEG$KO > 0.31315406,]
dim(data_DEG)



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

## 保存
write.csv(row.names(data_DEG)[loc],'/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_up/DEG_KO_WT_vs_KO_FAC_KO_FC_1.2.csv')


## 处理330个基因的富集结果
inputDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/metascape2/DEG_KO_WT_vs_KO_FAC_KO_FC_1.2/'
resDir <- inputDir
plotMetascapeBar(inputDir,resDir,goColor,keggColor)




## 寻找先下调后上调的基因
## 先取出部分基因
rm(list = ls())
k_40_res <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/kmeans_40/kmeans_result.csv')

## 读取Seurat对象
SeuratOb <- qread(paste0(RDSDir,'SeuratOb_UMAP_CR_zhuShi.qs'))


## 提取出cluster是 3 7 11 16 19 22 33 36 39的基因
k_40_res_use <- k_40_res %>%
  filter(cluster %in% c(3,7,11,16,19,22,33,36,39))

head(k_40_res_use)
use_genes <- k_40_res_use$rn
dim(k_40_res_use)

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



dim(SeuratOb_use)



## 跑差异基因（全部细胞类型，基因先下降后上调的基因）
umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_down/')
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
cellDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/ALL_KO_FAC_vs_KO_use_genes_down/')
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
DEG_KO_WT <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_down/ALL_markers.csv')
DEG_KO_FAC_KO <- fread('/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_KO_FAC_vs_KO_use_genes_down/ALL_markers.csv')

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
write.csv(row.names(data_DEG)[loc],'/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/ALL_use_genes_down/DEG_KO_WT_vs_KO_FAC_KO_FC_1.2.csv')


inputDir <- '/Users/lin/Desktop/backup/project/tfr1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/metascape2/down_DEG_KO_WT_vs_KO_FAC_KO_FC_1.2/'
resDir <- inputDir
plotMetascapeBar(inputDir,resDir,goColor,keggColor)
