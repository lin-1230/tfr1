######################################
# @date: 2025-04-09
# 这个脚本是为了重新对Tfr1课题的单细胞数据进行分群和注释
# 数据来源于刘宇进行去批次和整合后的单细胞数据
######################################

#### 在服务器环境中运行
# conda activate Tfr1
# R

## 设置当前工作目录
setwd("/media/ssd/sdb1/data/ljh/TFR1/")
# setwd("/Users/lin/Desktop/backup/project/tfr1/")

resultDir <- 'result/2025-04-09-zhuShi/'
dir.create(resultDir, recursive = TRUE)
RDSDir <- 'RDS/'
# RDSDir <- '/Volumes/julab/课题数据/TFR1/RDS'
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





## 读取刘宇处理的Seurat v5数据

SeuratOb <- qread('RDS/scRNA_merge_doublet.qs')
# SeuratOb <- readRDS(paste0(RDSDir,'/scRNA_merge_doublet.rds'))

colnames(SeuratOb@assays$RNA@layers$counts)

## 提取表达谱数据
expression_matrix <- LayerData(SeuratOb, layer = "counts")
meta.data <- SeuratOb@meta.data
harmony_embeddings <- Embeddings(SeuratOb, reduction = "harmony")

## 保存成qs文件
qsave(expression_matrix,paste0(RDSDir,'expression_matrix_raw.qs'))
qsave(meta.data,paste0(RDSDir,'meta.data_raw.qs'))
qsave(harmony_embeddings,paste0(RDSDir,'harmony_embeddings_raw.qs'))


## 创建Seurat v4版本的对象
SeuratOb2 <- CreateSeuratObject(counts = expression_matrix,meta.data=meta.data, min.cells = 3, min.features = 200)
head(SeuratOb2)

SeuratOb2[["harmony"]] <- CreateDimReducObject(embeddings = harmony_embeddings, key = "harmony_", assay = "RNA")
SeuratOb2[["harmony"]]

SeuratOb <- SeuratOb2
SeuratOb[['harmony']]


## 人的是MT-开头的基因
## 小鼠的是mt-开头的基因

## 计算线粒体基因的比例
SeuratOb[["percent.mt"]] <- PercentageFeatureSet(SeuratOb, pattern = "^mt-")



pdf(paste0(resultDir,'质控小提琴图.pdf'), width = pdf_width, height = pdf_height)
p1 <- VlnPlot(SeuratOb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3) &
 themeSet & 
 theme(axis.title = element_blank(),
        axis.text.x = element_blank())
print(p1)
dev.off()

pdf(paste0(resultDir,'质控小提琴图.pdf'), width = pdf_width, height = pdf_height)
p1 <- VlnPlot(SeuratOb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3) 
print(p1)
dev.off()

head(SeuratOb)


## 归一化
SeuratOb <- NormalizeData(SeuratOb)
## 寻找高变异基因，并绘制火山图
SeuratOb <- FindVariableFeatures(SeuratOb, selection.method = "vst", nfeatures = nfeatures)
top10 <- head(VariableFeatures(SeuratOb), 10)

pdf(paste0(resultDir,'高变异基因火山图.pdf'), width = 2 * pdf_width, height = pdf_height)
plot1 <- VariableFeaturePlot(SeuratOb) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# 调整图片中所有字体的大小
# 增加横坐标标签到图片的距离
print(plot1 + plot2 & themeSet & theme(axis.text.x = element_text(margin = margin(t = 15))))
dev.off()

SeuratOb <- ScaleData(SeuratOb, features = VariableFeatures(object = SeuratOb))

SeuratOb <- RunPCA(SeuratOb, features = VariableFeatures(object = SeuratOb),npcs=npcs)
print(SeuratOb[["pca"]], dims = 1:5, nfeatures = 20)

## 去掉SeuratOb对象中的scale.data 数据
## 减少内存占用，后续分析都不会用到scale.data
# SeuratOb@assays$RNA@scale.data <- matrix(nrow=0,ncol=0)

qsave(SeuratOb,paste0(RDSDir,'SeuratOb_PCA.qs'))
# SeuratOb <- qread(paste0(RDSDir,'SeuratOb_PCA.qs'))


## 寻找可能的最佳的pca维度
pdf(paste0(resultDir,'PCA基因贡献图.pdf'), width = pdf_width, height = pdf_height)
VizDimLoadings(SeuratOb, dims = 1:2, reduction = "pca") &
  themeSet 
dev.off()

pdf(paste0(resultDir,'维度选择图.pdf'), width = pdf_width, height = pdf_height)
ElbowPlot(SeuratOb) & 
  themeSet &
  theme(axis.text = element_text(size = axis_text_size))
dev.off()


## ----- PCA维度选择 -----
## 自动选择主成分数量算法：
## 1. 计算相邻主成分的标准差差值
## 2. 找到第一个差值小于阈值的维度
## 参数说明：dimShold=0.1 控制灵敏度
if(autoSelectDims){
  n <- length(SeuratOb@reductions[["pca"]]@stdev)
  a <- SeuratOb@reductions[["pca"]]@stdev[1:n-1]-SeuratOb@reductions[["pca"]]@stdev[2:n]
  selectPcaDim <- which(a < dimShold)[1]
}
selectPcaDim_raw <- selectPcaDim
print(paste0('选择的PCA纬度是：',selectPcaDim))


umapDir <- paste0(resultDir,'UMAP-harmony/')
for (resolution in c(0.1,0.3,0.5,0.7,0.9)){
  umapAndZhuShifun(SeuratOb = SeuratOb,
                   selectPcaDim = selectPcaDim,
                   resolution = resolution,
                   resultDir = umapDir,
                   reduction = 'harmony',
                   ifZhuShi = F,
                   ifplotsample = T)
}




umapDir <- paste0(resultDir,'UMAP-harmony/')
for (resolution in c(0.1,0.3,0.5,0.7,0.9)){
  umapAndZhuShifun(SeuratOb = SeuratOb,
                   selectPcaDim = 17,
                   resolution = resolution,
                   resultDir = umapDir,
                   reduction = 'harmony',
                   ifZhuShi = F,
                   ifplotsample = T)
}


################ 细胞注释
## 读取注释基因文件
markList <- fread('/media/ssd/sdb1/data/ljh/TFR1/genes/42003_2023_4936_MOESM3_ESM.csv')
colnames(markList)[1] <- 'LMPP_MPP'

head(markList)

metaData <- SeuratOb@meta.data
## 去掉cluster列
metaData <- metaData[,-which(colnames(metaData) %like% 'cluster')]
SeuratOb@meta.data <- metaData

SeuratOb <- umapAndZhuShifun(SeuratOb = SeuratOb,
                   selectPcaDim = 9,
                   resolution = 0.3,
                   resultDir = umapDir,
                   reduction = 'harmony',
                   ifZhuShi = T,
                   ifplotsample = T)




 
           
markList <- list(HSC=c("Hlf", "Mecom", "Ly6a","Procr"),
                MPP=c("Flt3", "Cd34", "Gcnt2", "Mki67"),
                CMP_MkP=c("Pf4", "Cd9", "Itga2b", "Angpt1", "Gata2"),
                MEP=c("Cd34", "Ly6a", "Fcgr3", "Gata1", "Gata1","Tfrc","Klf1"),
                GMP=c("Fcgr3", "F13a1", "Irf8", "Csf1r", "Mpo"),
                ProNeu=c("Cebpe", "Gfi1", "Prtn3", "Elane", "Mpo"),
                matureNue=c( "Ly6c2", "Cebpe", "Mpo"),
                DC=c("Irf8", "Klf4"),
                PC=c("Jchain", "Cd74"),
                NK_T=c('Lck', 'Ncr1'))

markList <- as.data.table(markList)


umapDir_DR <- paste0(resultDir,'UMAP-harmony_DR/')
SeuratOb <- umapAndZhuShifun(SeuratOb = SeuratOb,
                   selectPcaDim = 9,
                   resolution = 0.3,
                   resultDir = umapDir_DR,
                   reduction = 'harmony',
                   ifZhuShi = T,
                   ifplotsample = T)


## 注释结果
# 0 HSC
# 1 MPP
# 2 MPP
# 3 HSC
# 4 MPP
# 5 DC
# 6 PC
# 7 NK/T
# 8 MEP
# 9 

SeuratOb2 <- SeuratOb
 
head(SeuratOb)
Idents(SeuratOb)
table(Idents(SeuratOb))


new.cluster.ids <- c("HSC", "MPP", "MPP", "HSC", "MPP", "DC", 
                    "PC", "NK/T", "MEP", "unknown")

names(new.cluster.ids) <- levels(SeuratOb)
SeuratOb <- RenameIdents(SeuratOb, new.cluster.ids)


SeuratOb$group <- factor(SeuratOb$group,levels = c('WT1','WT2','WT_FAC','KO2','KO3','KO4','KO_FAC'))
levels(SeuratOb$group)




#### 绘制umap图
pdfDir <- paste0(umapDir_DR,'分群结果_dim-9_resolution-0.3/')
pdf(paste0(pdfDir,'聚类结果.pdf'),width=pdf_width,height=pdf_height)
DimPlot(SeuratOb,cols = col,pt.size = 1,label = T) & themeSet 
dev.off()
 


## 绘制每个分组的umap图
for (groupName in levels(SeuratOb$group)){
  plotGroupUmap(SeuratOb = SeuratOb,groupName = groupName,pdfDir = pdfDir,col = col,pt.size = 1,label = T)
}

table(SeuratOb$group)


SeuratOb$cellType <- Idents(SeuratOb)

### 细胞比例图
cell_data <- SeuratOb@meta.data %>%
  group_by(group, cellType) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>% # 计算比例
  ungroup()



pdf(paste0(pdfDir,'细胞比例图.pdf'),width = pdf_width,height=pdf_height)
ggplot(cell_data, aes(x = group, y = proportion, fill = factor(cellType))) +
  geom_bar(stat = "identity")  +
      labs(x = "Group", y = "Proportion", fill = "Cluster") +
      scale_fill_manual(values = col) & themeSet &
      theme(axis.text.x = element_text(angle = 0, hjust = 1))  # 使x轴标签倾斜，便于显示
dev.off()


### 细胞比例图（占全部细胞的占比）
cell_data <- SeuratOb@meta.data %>%
  group_by(group, cellType) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(proportion = count / sum(count))  # 计算比例
  

pdf(paste0(pdfDir,'全部细胞比例图.pdf'),width = pdf_width,height=pdf_height)
ggplot(cell_data, aes(x = group, y = proportion, fill = factor(cellType))) +
  geom_bar(stat = "identity")  +
      labs(x = "Group", y = "Proportion", fill = "Cluster") +
      scale_fill_manual(values = col) & themeSet &
      theme(axis.text.x = element_text(angle = 0, hjust = 1))  # 使x轴标签倾斜，便于显示
dev.off()






## 保存Seurat对象
qsave(SeuratOb,paste0(RDSDir,'SeuratOb_UMAP.qs'))
# SeuratOb <- qread(paste0(RDSDir,'SeuratOb_UMAP.qs'))

SeuratOb$cellType <- Idents(SeuratOb)

## 提取marker基因

markerGenes <- unique(unname(unlist(markList)))



levels(SeuratOb$cellType)


## 绘制细胞注释的marker基因点图（后续可能要继续优化）
DotPlot(SeuratOb, features = markerGenes, group.by = 'cellType',dot.scale=10) + ## dot.scale 定义点的大小
  coord_flip()+    ## 翻转坐标轴
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#5a67dc','#add8e6','#d6706d')) + 
  labs(x = "", y = "", title = "") +
  geom_tile(data = data.frame(cellType = levels(SeuratOb$cellType)), 
            aes(x = length(markerGenes) + 1, y = cellType, fill = cellType),  ## 修改色块的位置
            width = -1, height = 1.2, color = "white", inherit.aes = FALSE) + ## width是修改色块的高度，height是修改色块的宽度
  scale_fill_manual(values = col)  &  # 使用预定义的颜色映射
      themeSet + ## 加上全局的主题设置  + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())

  
# ggsave('result/聚类分群结果-env-seurat2/分群结果_dim-11_resolution-0.7/zhuShi/ALL/dot_markers_zhushi.pdf', width = 10, height = 6)
ggsave(paste0(pdfDir,'zhuShi/dot_markers_zhushi.pdf'), width = pdf_width, height = pdf_height)






## 使用CR文章中的marker基因进行注释（最后使用版本！）

markList=c('Ly6c2','Cebpe','Elane','Prtn3',
          'Vcam1','Gfi1','Mpo','Prss34','Fcer1a',
          'Ms4a2','F13a1','Klf4','Csf1r',
          'H2-Aa','Prdm1','Vpreb2','Lck',
          'Ncr1','Rag2','Dntt','Flt3','Mki67','Gata1',
          'Tspo2','Car1','Itga2b','Pf4','Ly6a','Cd48',
          'Cd34','Hlf','Clec1a','Procr','Sult1a1')

markList <- as.data.table(markList)


umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
SeuratOb <- umapAndZhuShifun(SeuratOb = SeuratOb,
                   selectPcaDim = 9,
                   resolution = 0.3,
                   resultDir = umapDir_CR,
                   reduction = 'harmony',
                   ifZhuShi = T,
                   ifplotsample = T)


## 按照师姐的要求，单独绘制了几个流式常用marker用来判断
markList = c('Cd34','Cd48','Flt3','Slamf1')
markList <- as.data.table(markList)

umapDir_hspc <- paste0(resultDir,'UMAP-harmony_hspc/')
SeuratOb <- umapAndZhuShifun(SeuratOb = SeuratOb,
                   selectPcaDim = 9,
                   resolution = 0.3,
                   resultDir = umapDir_hspc,
                   reduction = 'harmony',
                   ifZhuShi = T,
                   ifplotsample = T)


qsave(SeuratOb,paste0(RDSDir,'SeuratOb_UMAP_CR_zhuShi.qs'))
## 读取数据
# SeuratOb <- qread(paste0(RDSDir,'SeuratOb_UMAP_CR_zhuShi.qs'))


## 注释结果
# 0 HSC
# 1 cycling MPP2/3
# 2 MPP2/3
# 3 MPP1(ST)
# 4 MPP4
# 5 DC
# 6 PC
# 7 NK/T
# 8 Ery_P
# 9 B

## 注释结果
new.cluster.ids <- c("HSC", "cycling_MPP2/3", "MPP2/3", "MPP1", "MPP4", "DC",
                    "PC", "NK/T", "Ery_P", "B")

names(new.cluster.ids) <- levels(SeuratOb)
SeuratOb <- RenameIdents(SeuratOb, new.cluster.ids)
SeuratOb$cellType <- Idents(SeuratOb)

## 保存Seurat对象
qsave(SeuratOb,paste0(RDSDir,'SeuratOb_UMAP_CR_zhuShi.qs'))
## 读取数据
# SeuratOb <- qread(paste0(RDSDir,'SeuratOb_UMAP_CR_zhuShi.qs'))


pdfDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/')


DimPlot(SeuratOb,cols = col,pt.size = 1,label = T) & themeSet 
ggsave(paste0(pdfDir,'聚类结果.pdf'),width=pdf_width,height=pdf_height)


## TODO 绘制每个分组的umap图（讨论一下是不是要下采样，采样到多少？）
for (groupName in levels(SeuratOb$group)){
  plotGroupUmap(SeuratOb = SeuratOb,groupName = groupName,pdfDir = pdfDir,col = col,pt.size = 1,label = T)
}
table(SeuratOb$group)


SeuratOb$cellType <- Idents(SeuratOb)

### 细胞比例图
cell_data <- SeuratOb@meta.data %>%
  group_by(group, cellType) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>% # 计算比例
  ungroup()




ggplot(cell_data, aes(x = group, y = proportion, fill = factor(cellType))) +
  geom_bar(stat = "identity")  +
      labs(x = "Group", y = "Proportion", fill = "Cluster") +
      scale_fill_manual(values = col) & themeSet &
      theme(axis.text.x = element_text(angle = 0, hjust = 1))  # 使x轴标签倾斜，便于显示
ggsave(paste0(pdfDir,'细胞比例图.pdf'),width = pdf_width,height=pdf_height)


### 细胞比例图（占全部细胞的占比）
cell_data <- SeuratOb@meta.data %>%
  group_by(group, cellType) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(proportion = count / sum(count))  # 计算比例
  


ggplot(cell_data, aes(x = group, y = proportion, fill = factor(cellType))) +
  geom_bar(stat = "identity")  +
      labs(x = "Group", y = "Proportion", fill = "Cluster") +
      scale_fill_manual(values = col) & themeSet &
      theme(axis.text.x = element_text(angle = 0, hjust = 1))  # 使x轴标签倾斜，便于显示
ggsave(paste0(pdfDir,'全部细胞比例图.pdf'),width = pdf_width,height=pdf_height)




markerGenes <- unique(unname(unlist(markList)))

SeuratOb$cellType <- Idents(SeuratOb)
## 绘制细胞注释的marker基因点图（后续可能要继续优化）
DotPlot(SeuratOb, features = markerGenes, group.by = 'cellType',dot.scale=10) + ## dot.scale 定义点的大小
  coord_flip()+    ## 翻转坐标轴
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#5a67dc','#add8e6','#d6706d')) + 
  labs(x = "", y = "", title = "") +
  geom_tile(data = data.frame(cellType = levels(SeuratOb$cellType)), 
            aes(x = length(markerGenes) + 1, y = cellType, fill = cellType),  ## 修改色块的位置
            width = -1, height = 1.2, color = "white", inherit.aes = FALSE) + ## width是修改色块的高度，height是修改色块的宽度
  scale_fill_manual(values = col)  &  # 使用预定义的颜色映射
      themeSet + ## 加上全局的主题设置  + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())



  
# ggsave('result/聚类分群结果-env-seurat2/分群结果_dim-11_resolution-0.7/zhuShi/ALL/dot_markers_zhushi.pdf', width = 10, height = 6)
ggsave(paste0(pdfDir,'zhuShi/dot_markers_zhushi.pdf'), width = pdf_width, height = pdf_height)


## 绘制Tfr1的umap点图
'Tfrc' %in% rownames(SeuratOb)

umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
pdfDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/')

for (groupName in levels(SeuratOb$group)){
  ## 取出分组数据
  SeuratOb_use <- subset(SeuratOb, subset = group == groupName)

  pdf(paste0(pdfDir,groupName,'_Tfr1_umap.pdf'),width=pdf_width,height=pdf_height)
  print(FeaturePlot(SeuratOb_use,features = 'Tfrc',raster=FALSE) & themeSet)
  dev.off()
}


# ggsave(paste0(pdfDir,'Tfr1_umap.pdf'),width=pdf_width*3,height=pdf_height,limitsize=F)
    





## 取出HSC亚群进行重新分群，看看各个分组之间是否存在比例的差异
## 取出SHC细胞群
SeuratOb_HSC <- subset(SeuratOb, idents = c('HSC'))
table(SeuratOb_HSC$cellType)

## 进行重新分群

## 自动选取纬度
if(autoSelectDims){
  n <- length(SeuratOb_HSC@reductions[["pca"]]@stdev)
  a <- SeuratOb_HSC@reductions[["pca"]]@stdev[1:n-1]-SeuratOb_HSC@reductions[["pca"]]@stdev[2:n]
  selectPcaDim <- which(a < dimShold)[1]  # 找到第一个差值小于阈值的维度   
  
}
selectPcaDim_raw <- selectPcaDim
print(paste0('选择的PCA纬度是：',selectPcaDim))


umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR_HSC/')
umapRes <- list()
i <- 0
for (resolution in c(0.1,0.3,0.5,0.7,0.9)){
  i <- i + 1
  umapRes[i] <- umapAndZhuShifun(SeuratOb = SeuratOb_HSC,
                   selectPcaDim = selectPcaDim,
                   resolution = resolution,
                   resultDir = umapDir_CR,
                   reduction = 'harmony',
                   ifZhuShi = F,
                   ifplotsample = T)
}

## 绘制每个分组的umap图
resolutionV <- c(0.1,0.3,0.5,0.7,0.9)
for (i in 1:length(umapRes)){
  ## 分别合并WT重复和KO重复
  SeuratOb_use <- umapRes[[i]]
  SeuratOb_use$group2 <- as.character(SeuratOb_use$group)
  ## KO2、KO3、KO4 合并为KO
  SeuratOb_use$group2[SeuratOb_use$group2 %in% c('KO2','KO3','KO4')] <- 'KO'
  ## WT1、WT2合并为WT
  SeuratOb_use$group2[SeuratOb_use$group2 %in% c('WT1','WT2')] <- 'WT'

  

  ## 绘制堆砌条形图


  cell_data <- SeuratOb_use@meta.data %>%
    group_by(group2, seurat_clusters) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count)) %>% # 计算比例
    ungroup()

  # levels(cell_data$group2) <- c('WT','KO','WT_FAC','KO_FAC')
  cell_data$group2 <- factor(cell_data$group2,levels = c('WT','KO','WT_FAC','KO_FAC'))
  pdfDir <- paste0(umapDir_CR,'分群结果_dim-',selectPcaDim_raw,'_resolution-',resolutionV[i],'/')
  pdf(paste0(pdfDir,'聚类结果_堆砌条形图_合并重复.pdf'),width=8,height=8)
  p1 <- ggplot(cell_data, aes(x = group2, y = proportion, fill = factor(seurat_clusters))) +
      geom_bar(stat = "identity")  +
      labs(x = "Group", y = "Proportion", fill = "Cluster") +
      scale_fill_manual(values = col) & themeSet &
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 使x轴标签倾斜，便于显示
    print(p1)
  dev.off()
}






## Tfrc表达的比较
## 获取Tfr1表达值
SeuratOb$Tfrc <- FetchData(SeuratOb, vars = 'Tfrc')$Tfrc
head(SeuratOb)
## 按照分组计算Tfrc的均值
## 分组

library(data.table)
meta.data <- as.data.table(SeuratOb@meta.data)
meta.data <- meta.data[,.(Tfrc = mean(Tfrc)),by = .(group)]
table(SeuratOb$group)


umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
pdfDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/')
pdf(paste0(pdfDir,'Tfrc_boxplot.pdf'),width = pdf_width,height = pdf_height)
p1 <- ggplot(SeuratOb@meta.data, aes(x = group, y = Tfrc, fill = group))+
stat_summary(geom = "bar", fun = mean, position = position_dodge()) +  # 绘制均值条形图
stat_summary(geom = "errorbar",  # 添加误差棒
              fun.data = mean_se,  # 计算均值±标准误
              width = 0.2,         # 误差棒宽度
              ## 粗细
              linewidth = 1,
              position = position_dodge(0.9)) +  # 与条形图对齐
labs(y='Gene Expression')& themeSet
# geom_signif(comparisons = list(c("FuoAL", "NcoAL"),c("FuoDR", "NcoDR"),c("FuoDR", "FuoAL"),c("NcoDR", "NcoAL")),
#             map_signif_level = TRUE,
#             textsize = 7,
#             vjust = -0.5) 

print(p1)
dev.off()






## 每种细胞类型进行差异分析
SeuratOb <- qread(paste0(RDSDir,'SeuratOb_UMAP_CR_zhuShi.qs'))



## HSC


umapDir_CR <- paste0(resultDir,'UMAP-harmony_CR/')
hscDir <- paste0(umapDir_CR,'分群结果_dim-9_resolution-0.3/DEG/HSC/')
dir.create(hscDir, recursive = TRUE)
DEfunction(SeuratOb = SeuratOb,
           cellType = 'HSC',
           resDir = hscDir,
           group1 = 'KO',
           group2 = 'WT',
           logfc.threshold = 1,
           p_val.threshold = 0.05)




