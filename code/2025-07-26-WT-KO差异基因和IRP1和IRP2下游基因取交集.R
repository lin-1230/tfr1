## @date 2025-07-26
## @author lin
## @description WT\KO差异基因和IRP1和IRP2下游基因取交集


# 1. 进行差异分析
# 2. 读取差异分析结果
# 3. 读取IRP1和IRP2下游基因
# 4. 取交集
# 5. 输出结果到excel表
# 6. PPT整理

# conda activate RNAseq
# R

## 加载R包
library(Seurat)
library(ggplot2)
library(qs)
library(dplyr)
library(openxlsx)
library(readxl)
library(data.table)


setwd("/media/ssd/sdb1/data/ljh/TFR1/")


resultDir <- 'result/20250726/'
dir.create(resultDir, recursive = TRUE)
RDSDir <- 'RDS/'
dir.create(RDSDir, recursive = TRUE)


## 读取RDS数据
SeuratOb <- qread(paste0(RDSDir,'SeuratOb_UMAP_CR_zhuShi.qs'))
## check 一下数据
head(SeuratOb)

## 读取IRP1和IRP2下游基因
IRP <- read_excel('/media/ssd/sdb1/data/ljh/TFR1/genes/IRP.xlsx')
IRP <- IRP[-1,c(1,2)]
head(IRP)

## 读取差异基因
DEGdir='/media/ssd/sdb1/data/ljh/TFR1/result/2025-04-09-zhuShi/UMAP-harmony_CR/分群结果_dim-9_resolution-0.3/DEG/'

## 自定义读取并筛选差异基因的方程
# 定义读取差异基因的函数
readDEG <- function(DEGdir, DEclass, pvalue = 0.05, logfc = 1) {
  ## 参数解释
  # DEGdir: 差异基因目录
  # DEclass: 差异基因分类，分为ALL和具体的细胞类型
  # pvalue: 差异基因p值阈值
  # logfc: 差异基因logfc阈值

  # 读取差异基因文件
  DEG_all <- fread(paste0(DEGdir, DEclass, '/', DEclass, '_markers.csv'))
  print(head(DEG_all))
  
  # 筛选上调和下调的差异基因
  DEG_up <- DEG_all[avg_log2FC > logfc & p_val < pvalue,]
  DEG_down <- DEG_all[avg_log2FC < -logfc & p_val < pvalue,]
  
  cat('全部基因数量:', nrow(DEG_all), '\n')
  cat('上调基因数量:', nrow(DEG_up), '\n')
  cat('下调基因数量:', nrow(DEG_down), '\n')


  # 返回结果列表
  return(list(
    all = DEG_all,
    up = DEG_up,
    down = DEG_down
  ))
}

# 调用函数读取差异基因
# DEG_All_result <- readDEG(DEGdir, DEclass = 'ALL', pvalue = 0.05, logfc = 1)
## ALL
DEG_All_result_1 <- readDEG(DEGdir, DEclass = 'ALL', pvalue = 0.05, logfc = 1)
DEG_All_result_0p5 <- readDEG(DEGdir, DEclass = 'ALL', pvalue = 0.05, logfc = 0.5)
## HSC
DEG_HSC_result_1 <- readDEG(DEGdir, DEclass = 'HSC', pvalue = 0.05, logfc = 1)
DEG_HSC_result_0p5 <- readDEG(DEGdir, DEclass = 'HSC', pvalue = 0.05, logfc = 0.5)
## MPP1
DEG_MPP1_result_1 <- readDEG(DEGdir, DEclass = 'MPP1', pvalue = 0.05, logfc = 1)
DEG_MPP1_result_0p5 <- readDEG(DEGdir, DEclass = 'MPP1', pvalue = 0.05, logfc = 0.5)
## MPP2_3
DEG_MPP2_3_result_1 <- readDEG(DEGdir, DEclass = 'MPP2_3', pvalue = 0.05, logfc = 1)
DEG_MPP2_3_result_0p5 <- readDEG(DEGdir, DEclass = 'MPP2_3', pvalue = 0.05, logfc = 0.5)
## MPP4
DEG_MPP4_result_1 <- readDEG(DEGdir, DEclass = 'MPP4', pvalue = 0.05, logfc = 1)
DEG_MPP4_result_0p5 <- readDEG(DEGdir, DEclass = 'MPP4', pvalue = 0.05, logfc = 0.5)
## cycling_MPP2_3
DEG_cMPP2_3_result_1 <- readDEG(DEGdir, DEclass = 'cycling_MPP2_3', pvalue = 0.05, logfc = 1)
DEG_cMPP2_3_result_0p5 <- readDEG(DEGdir, DEclass = 'cycling_MPP2_3', pvalue = 0.05, logfc = 0.5)


## 提取差异基因
## ALL
## 创建两个list分别存储logfc>1和logfc>0.5的结果
DEG_list_1 <- list(
  ALL = list(
    up = DEG_All_result_1$up,
    down = DEG_All_result_1$down
  ),
  HSC = list(
    up = DEG_HSC_result_1$up,
    down = DEG_HSC_result_1$down
  ),
  MPP1 = list(
    up = DEG_MPP1_result_1$up,
    down = DEG_MPP1_result_1$down
  ),
  MPP2_3 = list(
    up = DEG_MPP2_3_result_1$up,
    down = DEG_MPP2_3_result_1$down
  ),
  MPP4 = list(
    up = DEG_MPP4_result_1$up,
    down = DEG_MPP4_result_1$down
  ),
  cMPP2_3 = list(
    up = DEG_cMPP2_3_result_1$up,
    down = DEG_cMPP2_3_result_1$down
  )
)

DEG_list_0p5 <- list(
  ALL = list(
    up = DEG_All_result_0p5$up,
    down = DEG_All_result_0p5$down
  ),
  HSC = list(
    up = DEG_HSC_result_0p5$up,
    down = DEG_HSC_result_0p5$down
  ),
  MPP1 = list(
    up = DEG_MPP1_result_0p5$up,
    down = DEG_MPP1_result_0p5$down
  ),
  MPP2_3 = list(
    up = DEG_MPP2_3_result_0p5$up,
    down = DEG_MPP2_3_result_0p5$down
  ),
  MPP4 = list(
    up = DEG_MPP4_result_0p5$up,
    down = DEG_MPP4_result_0p5$down
  ),
  cMPP2_3 = list(
    up = DEG_cMPP2_3_result_0p5$up,
    down = DEG_cMPP2_3_result_0p5$down
  )
)

## 自定义一个函数，分别与IRP1和IRP2 的下游基因取交集，并保存
DEGList <- DEG_list_1

interGeneSet <- function(DEGList , IRP, outDir){
  ## 参数解释
  ## DEGList: 差异分析结果列表
  ## IRP: IRP1和IRP2的下游基因列表
  ## outDir: 输出目录
  
  ## 遍历DEGList中每种细胞类型
  for (cellType in names(DEGList)){
    ## 分别取交集
    ## UP-IRP1
    res_up_IPR1 <- DEGList[[cellType]]$up[V1 %in% IRP$IRP1_targets,]
    ## UP-IRP2
    res_up_IPR2 <- DEGList[[cellType]]$up[V1 %in% IRP$IRP2_targets,]
    ## DOWN-IRP1
    res_down_IPR1 <- DEGList[[cellType]]$down[V1 %in% IRP$IRP1_targets,]
    ## DOWN-IRP2
    res_down_IPR2 <- DEGList[[cellType]]$down[V1 %in% IRP$IRP2_targets,]
    
    ## 按照avg_log2FC排序
    res_up_IPR1 <- res_up_IPR1[order(res_up_IPR1$avg_log2FC, decreasing = T),]
    res_up_IPR2 <- res_up_IPR2[order(res_up_IPR2$avg_log2FC, decreasing = T),]
    res_down_IPR1 <- res_down_IPR1[order(res_down_IPR1$avg_log2FC, decreasing = F),]
    res_down_IPR2 <- res_down_IPR2[order(res_down_IPR2$avg_log2FC, decreasing = F),]

    ## 输出交集后的基因数量
    cat(cellType, 'UP-IRP1:', nrow(res_up_IPR1), '\n')
    cat(cellType, 'UP-IRP2:', nrow(res_up_IPR2), '\n')
    cat(cellType, 'DOWN-IRP1:', nrow(res_down_IPR1), '\n')
    cat(cellType, 'DOWN-IRP2:', nrow(res_down_IPR2), '\n')

    ## 保存结果
    write.csv(res_up_IPR1, file = paste0(outDir, '/', cellType, '_up_IPR1.csv'), row.names = F)
    write.csv(res_up_IPR2, file = paste0(outDir, '/', cellType, '_up_IPR2.csv'), row.names = F)
    write.csv(res_down_IPR1, file = paste0(outDir, '/', cellType, '_down_IPR1.csv'), row.names = F)
    write.csv(res_down_IPR2, file = paste0(outDir, '/', cellType, '_down_IPR2.csv'), row.names = F)
    
  }
}

## FC > 1
dir.create(paste0(resultDir, '/DEG_FC_1'),recursive = T)
interGeneSet(DEG_list_1, IRP,paste0(resultDir, '/DEG_FC_1'))

## FC > 0.5
dir.create(paste0(resultDir, '/DEG_FC_0p5'),recursive = T)
interGeneSet(DEG_list_0p5, IRP,paste0(resultDir, '/DEG_FC_0p5'))