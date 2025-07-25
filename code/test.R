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
 

plotGroupUmap <- function(SeuratOb,groupName,pdfDir,col,pt.size = 1,label = T){
  ## 取出分组细胞
  SeuratOb_group <- subset(x = SeuratOb,subset = group == groupName)
  pdf(paste0(pdfDir,groupName,'_聚类结果.pdf'),width=pdf_width,height=pdf_height)
  p1 <- DimPlot(SeuratOb_group,cols = col,pt.size = pt.size,label = label) & themeSet
  print(p1)
  dev.off()
}

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



