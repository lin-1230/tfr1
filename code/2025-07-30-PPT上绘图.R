# conda activate RNAseq
# R
## 加载R包
library(Seurat)
library(ggplot2)
library(qs)
library(data.table)
library(dplyr)


setwd("/media/ssd/sdb1/data/ljh/TFR1/")
source('code/config_seurat.R')

RDSDir <- 'RDS/'
dir.create(RDSDir, recursive = TRUE)
resDir <- 'result/20250730/'
dir.create(resDir, recursive = TRUE)


## 读取SeuratOb
SeuratOb <- qread(paste0(RDSDir,'SeuratOb_20250728.qs'))
head(SeuratOb@meta.data)

cell_data <- SeuratOb@meta.data %>%
    group_by(group2,cellType) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count)) %>% # 计算比例
    ungroup()

# levels(cell_data$group2) <- c('WT','KO','WT_FAC','KO_FAC')
cell_data$group2 <- factor(cell_data$group2,levels = c('WT','KO','WT_FAC','KO_FAC'))
pdfDir <- resDir
pdf(paste0(pdfDir,'聚类结果_堆砌条形图_合并重复.pdf'),width=pdf_width,height=pdf_height)
p1 <- ggplot(cell_data, aes(x = group2, y = proportion, fill = factor(cellType))) +
    geom_bar(stat = "identity")  +
    labs(x = "Group", y = "Proportion", fill = "Cluster") +
    scale_fill_manual(values = col) & themeSet &
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 使x轴标签倾斜，便于显示
  print(p1)
dev.off()
