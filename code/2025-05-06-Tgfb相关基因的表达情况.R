##############################################
# @date: 2025-05-06
# 用来看Tgfb相关基因在HSC细胞中的表达情况

##############################################


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
