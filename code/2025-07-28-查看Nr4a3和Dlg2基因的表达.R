######################################
# @date: 2025-07-28
# 这个脚本是为了对Tfr1的单细胞数据中Nr4a3和Dlg2基因的表达
######################################

#### 在服务器环境中运行
# conda activate RNAseq
# R

## 设置当前工作目录
setwd("/media/ssd/sdb1/data/ljh/TFR1/")

RDSDir <- 'RDS/'
dir.create(RDSDir, recursive = TRUE)


## 加载R包
library(Seurat)
# library(data.table)
# library(stringr)
# library(Matrix)
# library(ggplot2)
library(qs)
# library(dplyr)


## 引入配置文件




## 读取数据


SeuratOb$group2 <- as.character(SeuratOb$group)

## KO2、KO3、KO4 合并为KO
SeuratOb$group2[SeuratOb$group %in% c('KO2','KO3','KO4')] <- 'KO'
## WT1、WT2合并为WT
SeuratOb$group2[SeuratOb$group %in% c('WT1','WT2')] <- 'WT'
  
SeuratOb$group2 <- factor(SeuratOb$group2, levels = c('WT','WT_FAC','KO','KO_FAC'))

qsave(SeuratOb, file = paste0(RDSDir,'SeuratOb_20250728.qs'))


# 计算每个基因在每个组中的平均表达值
avg_expr <- AverageExpression(SeuratOb, 
                            features = c('Nr4a1','Nr4a3','Dlg2'),
                            group.by = "group2",
                            slot = "data")
# 将矩阵数据转换为热图格式
plot_matrix <- as.matrix(avg_expr$RNA)

# 对每一行进行scale标准化处理
plot_matrix_scaled <- t(scale(t(plot_matrix)))



# 使用pheatmap绘制热图
library(pheatmap)

pdf(file=paste0('result/20250728/','Nr4a1_Nr4a3_Dlg2_avgExpr.pdf'),height=2.5,width=6)
pheatmap(plot_matrix_scaled,
         display_numbers = TRUE,  # 显示具体数值
         cluster_rows = FALSE,    # 不对行进行聚类
         cluster_cols = FALSE,    # 不对列进行聚类
         main = "mean Exp of Nr4a1, Nr4a3 and Dlg2",
         fontsize = 10,          # 设置字体大小
         angle_col = 45,         # 设置列名角度
         number_format = "%.3f")  # 设置数值格式为3位小数
dev.off()



# pdf(file=paste0('result/20250728/','Tfrc_umap.pdf'),height=10,width=30)
# FeaturePlot(SeuratOb, features = c('Tfrc'), pt.size = 1,label = TRUE,split.by = 'group2')
# dev.off()

# AverageExpression(SeuratOb, 
#                             features = c('Tfrc'),
#                             group.by = "group2",
#                             slot = "data")


# AverageExpression(SeuratOb, 
#                             features = c('Tfrc'),
#                             group.by = "group2",
#                             slot = "counts")





## 查看Bnip3、Ppargc1a、Il10、Socs3、Ccl3、Cdkn1a 的表达情况



plotGeneExpressionHeatmap <- function(seurat_obj, 
                                    features,
                                    group_by,
                                    output_dir,
                                    output_filename,
                                    main,
                                    height = 2.5,
                                    width = 6) {
    

    ## 参数解释
    # seurat_obj: seurat对象
    # features: 要查看的基因列表
    # group_by: 分组变量
    # output_dir: 输出目录
    # output_filename: 输出文件名
    # main: 热图标题
    # height: 热图高度
    # width: 热图宽度


    # 计算平均表达值
    avg_expr <- AverageExpression(seurat_obj, 
                                features = features,
                                group.by = group_by,
                                slot = "data")
    
    # 将矩阵数据转换为热图格式
    plot_matrix <- as.matrix(avg_expr$RNA)
    
    # 对每一行进行scale标准化处理
    plot_matrix_scaled <- t(scale(t(plot_matrix)))
    
    # 使用pheatmap绘制热图
    library(pheatmap)
    
    pdf(file = file.path(output_dir, output_filename), height = height, width = width)
    pheatmap(plot_matrix_scaled,
             display_numbers = TRUE,  # 显示具体数值
             cluster_rows = FALSE,    # 不对行进行聚类
             cluster_cols = FALSE,    # 不对列进行聚类
             main = main,
             fontsize = 10,          # 设置字体大小
             angle_col = 45,         # 设置列名角度
             number_format = "%.3f")  # 设置数值格式为3位小数
    dev.off()
}

# 调用函数

plotGeneExpressionHeatmap(SeuratOb,
                         features = c('Bnip3','Ppargc1a','Il10','Socs3','Ccl3','Cdkn1a'),
                         group_by = "group2",
                         output_dir = "result/20250728",
                         output_filename = 'Bnip3_Ppargc1a_Il10_Socs3_Ccl3_Cdkn1a_avgExpr.pdf',
                         main = "mean Exp of Bnip3, Ppargc1a, Il10, Socs3, Ccl3, Cdkn1a",
                         height = 3,
                         width = 6)
