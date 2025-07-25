#####################################
## 来自CHOA课题
## date ： 2025-04-09
#####################################s



col = c('#5668E4','#FFA9AF','#E66868','#969696','#4497E5',
        '#AFE4EA','#F8CD42','#C486CF','#CC19EF','#429342',
        '#5DFCD5','#FFECB3FF','#FCF78A','#F59D30',
        '#B59C93','#C6312C','#8994AC','#254365','#BAC0F2',
        '#E7E7E7','#A2CCEA','#835C27','#433825','#CFA670')

## 自定义一个用来分群和聚类的函数
## @ date 2025/04/09 增加了reduction参数，默认使用pca
umapAndZhuShifun <- function(SeuratOb,resultDir,selectPcaDim,resolution,ifZhuShi,ifplotsample,pt.size=0.7,pdf_width=10,Neighbors = 30,pdf_height=10,reduction='pca'){
  print('--------------------------------------------------------------------')
  print('正在绘制：')
  print(paste0('dim:',selectPcaDim))
  print(paste0('resolution:',resolution))
  
  ## 根据维度和分辨率新建结果文件夹
  umap_dir <- paste0(resultDir,'分群结果_dim-',selectPcaDim,'_resolution-',resolution,'/')
  dir.create(umap_dir,recursive = T)
  
  SeuratOb <- FindNeighbors(SeuratOb, dims = 1:selectPcaDim,reduction = reduction)
  SeuratOb <- FindClusters(SeuratOb, resolution = resolution)
  
  
  SeuratOb <- RunUMAP(SeuratOb, dims = 1:selectPcaDim,reduction=reduction,n.neighbors = Neighbors)
  
  #sampleName <- colnames(SeuratOb)
  #sampleName <- str_split(sampleName,'_',simplify = T)[,1]
  
  
  
  ## 提取样本编号
  #SeuratOb$sampleName <- sampleName
  
  n <- length(levels(SeuratOb$seurat_clusters))
  plotcol <- col[1:n]
  


  ## 绘制UMAP图
  pdf(paste0(umap_dir,'聚类结果UMAP.pdf'),width = pdf_width,height = pdf_height)
  p1 <- DimPlot(SeuratOb, reduction = "umap",cols = col[1:n],label = T,pt.size = pt.size) & 
    themeSet
  print(p1)
  dev.off()
  


  ## @date 2024/09/08 新增绘制，每个样本分别的点图
  ## 这样可以更清楚地分辨每个样本主要在哪些区域富集
  ## 是否有批次效应
  ## 对样本名称进行去重
  if (ifplotsample){
    uni_group <- unique(SeuratOb$group)
    
    for (i in 1:length(uni_group)){
      s <- uni_group[i]
      # 提取相应的样本
      seurat_sub <- subset(SeuratOb, subset = group == s )
      pdf(paste0(umap_dir,'聚类结果UMAP_',s,'.pdf'),width = pdf_width,height = pdf_height)
      p1 <- DimPlot(seurat_sub, reduction = "umap",cols = col[1:n],label = T,pt.size = pt.size) &
        themeSet
      #p2 <- DimPlot(seurat_sub,reduction = "umap",group.by = 'group',cols = col[i])
      print(p1&themeSet)
      dev.off()
    

      write.csv(table(Idents(seurat_sub)),paste0(umap_dir,s,'_样本分群统计表.csv'))
    }

    ########### 绘制每个分组的细胞比例图#############
    library(dplyr)
    cell_data <- SeuratOb@meta.data %>%
      group_by(group, seurat_clusters) %>%
      summarise(count = n()) %>%
      mutate(proportion = count / sum(count)) %>% # 计算比例
      ungroup()

    print(cell_data)

    pdf(paste0(umap_dir,'聚类结果_堆砌条形图.pdf'), width = 8, height = 8)
    p1 <- ggplot(cell_data, aes(x = group, y = proportion, fill = factor(seurat_clusters))) +
      geom_bar(stat = "identity")  +
      labs(x = "Group", y = "Proportion", fill = "Cluster") +
      scale_fill_manual(values = col) & themeSet &
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 使x轴标签倾斜，便于显示
    print(p1)
    dev.off()

    ## 绘制分组umap图
    pdf(paste0(umap_dir,'聚类结果_分组UMAP图.pdf'),width = pdf_width,height = pdf_height)
    p1 <- DimPlot(SeuratOb, reduction = "umap",label=T,split.by = 'group',cols = col[1:length(levels(SeuratOb$seurat_clusters))])
    print(p1 & themeSet)
    dev.off()

    
    return(SeuratOb)

  }

    

  if (ifZhuShi){
    ## 设置注释结果存放地址
    zhuShi_dir <- paste0(umap_dir,'zhuShi/')
    dir.create(zhuShi_dir,recursive = T)
    
    zhuShi(SeuratOb,                  ## 单细胞数据
           outPutDir=zhuShi_dir,             ## 结果路径
           markList,                        ## marker的list
           ifFindMark=F,                    ## 是否运行findAllMarker
           ifSingleR=F,                   ## 是否运行singleR
           pdf_width=pdf_width,             ## 图片宽度
           pdf_height=pdf_height)           ## 图片高度            
  }
  
  return(SeuratOb)
  
}




############################################
# @date:2025/01/15
# 
# 该函数是单细胞数据中，用来 做自定义细胞注释出图的函数
# 2024/01/15 针对selp—hsc课题改良
#
############################################

zhuShi <- function(seuratOB,
                   outPutDir,
                   markList,
                   ifFindMark=F,
                   ifSingleR=F,
                   minPct=0.1,
                   logfcThreshold=0.25,
                   returnThresh=0.05,
                   topN=20,
                   pdf_width=10,
                   pdf_height=10){
  
  ## seuratOB,seurat对象
  ## minPct，Seurat包中findAllMarker函数中的
  ## min.pct参数，默认0.25
  ## 用来制定至少在多少比例中的细胞中表达的基因才会
  ## 被选为marker基因
  
  ## markList,data.table,包含了marker基因的data.table
  ## 每一列代表一个细胞类型，每一行代表一个marker基因
  ## 或者是一个list，每个list的元素是一个data.table
  ## 或者是一个向量
  
  ## outPutDir,character,所有结果输出的位置
  ## ifFindMark、ifSingleR，是否开启findMark、singleR流程
  ## 默认开启
  
  ## logfcThreshold、returnThresh,numeric,Seurat包中findAllMarker函数中的参数
  ## 分别用来限制logfc和p的阈值
  
  ## topN,每个Cluster提取的top Marker基因数量
  
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  
  
  dir.create(outPutDir,recursive = T)
  
  
  ##########寻找差异表达基因
  if (ifFindMark){
    markers <- FindAllMarkers(seuratOB, only.pos = TRUE, min.pct = minPct, logfc.threshold = logfcThreshold,return.thresh = returnThresh)
    #每个cluster的top基因
    top <- markers %>%
      group_by(cluster) %>%
      top_n(n = topN, wt = avg_log2FC)
    
    
    pdf(file=paste0(outPutDir,'heatMap.pdf'),width=pdf_width,height=pdf_height)
    print(DoHeatmap(seuratOB,                             #seurat对象(数据)
                    features = top$gene,            #绘制的基因
                    label = F,                        #图中是否添加label信息
                    slot = "scale.data",              #绘图使用的表达矩阵
                    group.by = 'seurat_clusters',     #分组名称
                    group.bar = T) & themeSet)                 #是否展示bar
    dev.off()
    write.table(top,paste0(outPutDir,'heatMap.xls'),sep='\t',quote = T,row.names=F)
    
  }
  
  
  
  #####################singleR注释
  if (ifSingleR){
    library(SingleR)
    library(celldex)
    library(Seurat)
    library(pheatmap)
    
    #seuratOB的meta文件，包含了seurat的聚类结果
    meta<-seuratOB@meta.data 
    
    #载入参考数据集
    ref_use <- HumanPrimaryCellAtlasData()
    #SingleR注释
    pred <- SingleR(test=as.matrix(seuratOB@assays$RNA@count),       #输入表达矩阵
                    ref=ref_use,                                #参考数据
                    labels=ref_use$label.fine,                  #标签列
                    clusters = seuratOB$seurat_clusters,         #细胞聚类信息（只有method = "cluster"才使用）
                    method = "cluster")                         #注释类型，按照cluster还是按照cell进行注释
    
    #保存SingleR注释结果
    write.table(pred,paste0(outPutDir,'singleR.xls'),sep='\t',quote = F)
  }
  
  
  #####################根据marker注释
  
  ## @date 2024/09/10
  ## 加入所有细胞marker的气泡图，更加直观的观看
  ## marker在不同群细胞中的表达变化（还没测试）
  
  if (class(markList)[1] != 'data.table'){
    message('marker必须输入一个data.table对象！')
    return()
  }
  
  cellTypeNameList <- colnames(markList)
  for (i in 1:ncol(markList)){
    marker <- markList[,i,with = F]
    marker <- unique(marker)
    cellTypeName <- cellTypeNameList[i]
    ## 创建结果文件夹
    dir.create(paste0(outPutDir,'/',cellTypeName,'/'),recursive = T)
    
    DotPlot(seuratOB, features = marker, cluster.idents=T)+coord_flip()+
      theme_bw()+
      theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
      labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
      #scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
      scale_color_gradientn(values = seq(0,1,0.2),colours = c('white','#330066')) &
      themeSet  ## 加上全局的主题设置
    ggsave(paste0(outPutDir,'/',cellTypeName,'/dot_markers_',cellTypeName,'.pdf'),width = pdf_width,height = pdf_height)
    
  }
  
  
  ## 画一个总的点图
  markVer <- unique(as.character(unlist(markList)))
  dir.create(paste0(outPutDir,'/ALL/'),recursive = T)
  
  DotPlot(seuratOB, features = markVer, cluster.idents=T)+coord_flip()+
    theme_bw()+
    theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
    #scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('white','#330066')) & 
    themeSet  ## 加上全局的主题设置
  ggsave(paste0(outPutDir,'/ALL/dot_markers_ALL.pdf'),width = pdf_width,height = pdf_height)
  
  ## 点图
  dir.create(paste0(outPutDir,'/FeaturePlot/'),recursive = T)
  for (i in 1:length(markVer)){
    marker <- markVer[i]
    FeaturePlot(seuratOB,features = marker,raster=FALSE) & themeSet
    ggsave(paste0(outPutDir,'/FeaturePlot/',marker,'.pdf'),width = pdf_width,height = pdf_height)
    
  }
} 








## 自定义函数用于给HSC两群亚群作差异分析
## 并返回差异分析的结果，以及Selp在其中的排名情况
findSelp <- function(SeuratOb_HSC,test.use = 'wilcox',ifreturnMarkerDF=F,min.pct = 0.05,logfc.threshold = 0){
  print(paste0('现在正在使用的差异方法是：',test.use))
  ## 查看两群细胞的marker基因
  markers <- FindAllMarkers(SeuratOb_HSC, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold,return.thresh = 0.05,test.use = test.use)
  markers <- as.data.table(markers)
  print('Selp 在差异结果中的情况：')
  print(markers[gene=='Selp'])
  markers_1 <- markers[cluster=='1']
  ## 按照p_val排序，并给出排名
  if (test.use == 'roc'){
    markers_1 <- setorder(markers_1,myAUC)
  }else{
    markers_1 <- setorder(markers_1,p_val)
  }
  markers_1$rank <- row.names(markers_1)
  ## 查看Selp的排名以及p值
  print('Selp在其中的排名情况：')
  print(markers_1[gene=='Selp',rank])
  print('c1总共的marker基因数目是')
  print(nrow(markers_1))
  print('----------------------------------------------')
  if (ifreturnMarkerDF){
    return(markers)
  }else{
    return(markers_1$gene) 
  }
}


## 自定义一个用来给HSC进行MK基因打分的函数
mkScore <- function(SeuratOb_HSC,mkGenes,geneListName,pdfDir){
  
  library(viridis)
  library(ggplot2)
  library(ggsignif)
  
  
  ## 对HSC亚群进行打分
  SeuratOb_HSC_mkScore <- AddModuleScore(object = SeuratOb_HSC,features = list(mkGenes))
  
  ## 先对打分进行归一化处理
  SeuratOb_HSC_mkScore$Cluster1 <- scale(SeuratOb_HSC_mkScore$Cluster1)
  
  
  #pdf(paste0('result/seuratPDF/分群结果_dim-11_resolution-0.6/','HSC_打分_MKP.pdf'),height = 10)
  #p1 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Cluster1',shape.by = 'group',raster=FALSE,pt.size = 5,cols = c('#E7E7E7','#FFA9AF','#E66868','#F59D30'))
  #p2 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Cluster1',shape.by = 'seurat_clusters',raster=FALSE,pt.size = 5,cols = c('#E7E7E7','#FFA9AF','#E66868','#F59D30'))
  #print(p1+p2)
  #dev.off()
  
  pdf(paste0(pdfDir,'HSC_打分_',geneListName,'.pdf'),width = 15)
  p1 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Cluster1',split.by = 'group',shape.by = 'seurat_clusters',raster=FALSE,pt.size = 5,cols = viridis(6),label = T,
                    min.cutoff = -3,max.cutoff = 3)+scale_alpha(range = c(0.6, 1))
  
  p1 <- p1+ theme(legend.position = 'right')
  print(p1)
  dev.off()
  
  
  ## 画一个箱型图看看
  pdf(paste0(pdfDir,'HSC_打分_',geneListName,'_box.pdf'),width = 15)
  
  p1 <- ggplot(SeuratOb_HSC_mkScore@meta.data, aes(x = seurat_clusters, y = Cluster1, fill = seurat_clusters))+
    geom_boxplot()+
    labs(y='score')+
    geom_signif(comparisons = list(c('1','0')),  ##添加显著性
                map_signif_level = TRUE, 
                textsize = 4, 
                vjust = -0.5)
  print(p1)
  dev.off()
  
  
  pdf(paste0(pdfDir,'HSC_打分_Selp.pdf'),,width = 10)
  p1 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Selp',split.by = 'group',shape.by = 'seurat_clusters',raster=FALSE,pt.size = 5,cols = viridis(6),label = T,
                    min.cutoff = -3,max.cutoff = 3)+scale_alpha(range = c(0.6, 1))
  
  p1 <- p1+ theme(legend.position = 'right')
  print(p1)
  dev.off()
  
  pdf(paste0(pdfDir,'HSC_打分_Selp_vio.pdf'),,width = 10)
  p1 <- VlnPlot(SeuratOb_HSC_mkScore, features = "Selp", group.by = "seurat_clusters",pt.size = 1)
  print(p1)
  dev.off()
}



# 自动调整 y_position 的函数
generate_y_positions <- function(comparisons, data, offset = 1, max_offset = 10) {
  # 获取数据的最大值
  max_y <- max(data$Cluster1, na.rm = TRUE)
  
  # 计算每个显著性标记的 y_position，避免重叠
  y_positions <- seq(max_y + offset, max_y + offset * length(comparisons), by = offset)
  
  # 如果y_position超过最大偏移量，进行限制
  y_positions <- pmin(y_positions, max_y + max_offset)
  
  return(y_positions)
}



## 自定义用于多个群进行打分的函数
mkScoreMuti <- function(SeuratOb_HSC,mkGenes,geneListName,pdfDir,pt.size=2,boxOrder=NULL,boxPdfWidth=15){
  
  library(Seurat)
  library(viridis)
  library(ggplot2)
  library(ggsignif)
  
  
  ## 对HSC亚群进行打分
  SeuratOb_HSC_mkScore <- AddModuleScore(object = SeuratOb_HSC,features = list(mkGenes))
  
  ## 保存一份打分出来，有时候分析要用
  score <- SeuratOb_HSC_mkScore@meta.data
  
  
  ## 先对打分进行归一化处理
  #SeuratOb_HSC_mkScore$Cluster1 <- log10(SeuratOb_HSC_mkScore$Cluster1 + 10)
  SeuratOb_HSC_mkScore$Cluster1 <- scale(SeuratOb_HSC_mkScore$Cluster1)
  
  
  #pdf(paste0('result/seuratPDF/分群结果_dim-11_resolution-0.6/','HSC_打分_MKP.pdf'),height = 10)
  #p1 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Cluster1',shape.by = 'group',raster=FALSE,pt.size = 5,cols = c('#E7E7E7','#FFA9AF','#E66868','#F59D30'))
  #p2 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Cluster1',shape.by = 'seurat_clusters',raster=FALSE,pt.size = 5,cols = c('#E7E7E7','#FFA9AF','#E66868','#F59D30'))
  #print(p1+p2)
  #dev.off()
  
  pdf(paste0(pdfDir,'HSC_打分_',geneListName,'.pdf'),width = 15)
  p1 <- FeaturePlot(SeuratOb_HSC_mkScore,
                    features = 'Cluster1',
                    split.by = 'group',
                    raster=FALSE,
                    pt.size = pt.size,
                    cols = inferno(6),
                    label = T)+scale_alpha(range = c(0.6, 1))
  
  
  p1 <- p1+ theme(legend.position = 'right')
  print(p1)
  dev.off()
  
  
  
  pdf(paste0(pdfDir,'HSC_打分_',geneListName,'_whole.pdf'))
  p1 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Cluster1',
                    raster=FALSE,
                    pt.size = pt.size,
                    cols = inferno(6),
                    label = T)+
    scale_alpha(range = c(0.6, 1))
  
  
  p1 <- p1+ theme(legend.position = 'right')
  print(p1)
  dev.off()
  
  
  ## 画一个箱型图看看
  pdf(paste0(pdfDir,'HSC_打分_',geneListName,'_box.pdf'),width = boxPdfWidth)
  
  n <- length(unique(Idents(SeuratOb_HSC_mkScore)))
  #comparisons <- combn(0:n, 2,simplify = F)
  
  # 生成自动 y_position
  #_positions <- generate_y_positions(comparisons, SeuratOb_HSC_mkScore@meta.data, offset = 1, max_offset = 50)
  #print(y_positions)
  
  df <- SeuratOb_HSC_mkScore@meta.data
  
  ## 如果指定了顺序，那就按照指定的顺序来
  if(!is.null(boxOrder)){
    df$seurat_clusters  <- factor(df$seurat_clusters,levels = boxOrder)
  }else{
    # 按照 Cluster1 的值对 seurat_clusters 排序
    df$seurat_clusters <- reorder(df$seurat_clusters, df$Cluster1,decreasing=T, FUN = median)
  }
  
  ## 将seurat_clusters顺序返回到Seurat对象中，后面画Selp的小提琴图的时候可以用
  SeuratOb_HSC_mkScore$seurat_clusters <- df$seurat_clusters
  
  p1 <- ggplot(df, aes(x = seurat_clusters, y = Cluster1, fill = seurat_clusters))+
    geom_boxplot()+
    scale_fill_manual(values = col[1:n])+
    labs(y='score')+
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.25),  # 设置面板背景
          plot.title = element_text(hjust = 0.5),  # 设置标题水平居中对齐
          axis.line = element_line(colour = "black", linewidth = 0.25),  # 设置坐标轴线样式
          axis.title = element_text(size = 13, face = "plain", color = "black"),  # 设置坐标轴标题样式
          axis.text = element_text(size = 15, face = "plain", color = "black"),  # 设置坐标轴标签样式
          legend.text = element_text(face = "plain", colour = "black", size = 15),  # 设置图例文本样式
          legend.title = element_text(face = "plain", colour = "black", size = 15),  # 设置图例标题样式
          panel.grid.major = element_blank(),  # 隐藏主要网格线
          panel.grid.minor = element_blank(),  # 隐藏次要网格线
          axis.text.x = element_text(angle = 0, hjust = 0.5))  # 居中对齐
  
  #geom_signif(comparisons = comparisons,  ##添加显著性
  #            map_signif_level = TRUE, 
  #            textsize = 4,
  #            vjust = -0.5,
  #            y_position = y_positions)
  print(p1)
  dev.off()
  
  
  pdf(paste0(pdfDir,'HSC_打分_Selp.pdf'),width = 10)
  p1 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Selp',split.by = 'group',raster=FALSE,pt.size = pt.size,cols = c("#e4e3e3", "#b80000"),label = T,
                    min.cutoff = 0,max.cutoff = 1.5)+scale_alpha(range = c(0.6, 1))
  
  p1 <- p1+ theme(legend.position = 'right')
  print(p1)
  dev.off()
  
  pdf(paste0(pdfDir,'HSC_打分_Selp_whole.pdf'),width = 10)
  p1 <- FeaturePlot(SeuratOb_HSC_mkScore,features = 'Selp',raster=FALSE,pt.size = pt.size,cols = c("#e4e3e3", "#b80000"),label = T,
                    min.cutoff = 0,max.cutoff = 1.5)+scale_alpha(range = c(0.6, 1))
  
  p1 <- p1+ theme(legend.position = 'right')
  print(p1)
  dev.off()
  
  
  ## 绘制小提琴图
  pdf(paste0(pdfDir,'HSC_打分_Selp_vio.pdf'),,width = 10)
  p1 <- VlnPlot(SeuratOb_HSC_mkScore, features = "Selp", group.by = "seurat_clusters",pt.size = 1,cols = col[1:n])
  print(p1)
  dev.off()
  
  ## 用ggplot2 绘制 Selp图 
  df$rna_Selp <- FetchData(SeuratOb_HSC_mkScore,vars = 'rna_Selp')$rna_Selp
  
  
  pdf(paste0(pdfDir,'HSC_打分_Selp_box.pdf'),width = boxPdfWidth)
  p_selp_box <- ggplot(df, aes(x = seurat_clusters, y = rna_Selp, fill = seurat_clusters))+
    geom_boxplot()+
    scale_fill_manual(values = col[1:n])+
    labs(y='score')+
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.25),  # 设置面板背景
          plot.title = element_text(hjust = 0.5),  # 设置标题水平居中对齐
          axis.line = element_line(colour = "black", linewidth = 0.25),  # 设置坐标轴线样式
          axis.title = element_text(size = 13, face = "plain", color = "black"),  # 设置坐标轴标题样式
          axis.text = element_text(size = 15, face = "plain", color = "black"),  # 设置坐标轴标签样式
          legend.text = element_text(face = "plain", colour = "black", size = 15),  # 设置图例文本样式
          legend.title = element_text(face = "plain", colour = "black", size = 15),  # 设置图例标题样式
          panel.grid.major = element_blank(),  # 隐藏主要网格线
          panel.grid.minor = element_blank(),  # 隐藏次要网格线
          axis.text.x = element_text(angle = 0, hjust = 0.5))  # 居中对齐
  print(p_selp_box)
  dev.off()
  
  
  return(score)
  
}



## 自定义生成雷达图数据的函数
dealRadarData <- function(scoreM,pathwayName){
  library(data.table)
  
  meanScore <- list()
  for ( i in 1:length(scoreM)){
    scoreDF <- as.data.table(scoreM[[i]])
    
    scoreDF[, (pathwayName[i]) := mean(Cluster1), by = seurat_clusters]
    # 提取所需列（假设只需要 seurat_clusters 和新均值列）
    meanScore[[i]] <- scoreDF[, .SD, .SDcols = c("seurat_clusters", pathwayName[i]) ]
    meanScore[[i]] <- unique.data.frame(meanScore[[i]])
    
  }
  merged_dt <- Reduce(function(dt1, dt2) {
    merge(dt1, dt2, by = names(dt1)[1], all = TRUE)
  }, meanScore)
  
  return(merged_dt)
  
}


metascape2barplot <- function(upResFile,downResFile,pdfDir,colors=c('#4497E5','#E66868')){
  library(readxl)
  library(stringr)
  library(ggplot2)
  
  ## 处理上游结果
  #upResFile <- '/Users/lin/Desktop/backup/project/shaoTong/OA/result/差异/差异_SW-M1-HTvsWT/up_ALL/metascape_result.xlsx'
  up <- read_excel(upResFile,sheet = 2)
  loc <- str_split(up$GroupID,'_',simplify = T)[,2]=='Summary'
  up <- up[loc,]
  up$color <- rep('Up',nrow(up))
  up$LogP <- -1*up$LogP
  
  ## 处理下游结果
  #downResFile <- '/Users/lin/Desktop/backup/project/shaoTong/OA/result/差异/差异_SW-M1-HTvsWT/down_ALL/metascape_result.xlsx'
  down <- read_excel(downResFile,sheet = 2)
  loc <- str_split(down$GroupID,'_',simplify = T)[,2]=='Summary'
  down <- down[loc,]
  down$color <- rep('Stable',nrow(down))
  
  ## 合并数据
  barData <- rbind(down,up)
  
  
  ## 如果有重复
  if(sum(duplicated(barData$Description))){
    loc<-duplicated(barData$Description)
    barData$Description[loc] <- paste0(barData$Description[loc],'_2')
  }
  barData$Description <- factor(barData$Description,levels = barData$Description)
  
  
  pdf(pdfDir,width = 15)
  p1 <- ggplot(barData,aes(x = Description,y = LogP))+
    geom_col(aes(fill=color))+
    theme(axis.text.x = element_text(angle = 90))+ ## 横坐标标签调整角度
    ylab('log(p value)')+
    coord_flip() + 
    scale_fill_manual(values = colors) +
    theme(#legend.position = 'none',
      panel.background = element_rect(fill = "white", colour = "black", size = 0.25),  # 设置面板背景
      plot.title = element_text(hjust = 0.5),  # 设置标题水平居中对齐
      axis.line = element_line(colour = "black", linewidth = 0.25),  # 设置坐标轴线样式
      axis.title = element_text(size = 13, face = "plain", color = "black"),  # 设置坐标轴标题样式
      axis.text = element_text(size = 13, face = "plain", color = "black"),  # 设置坐标轴标签样式
      legend.text = element_text(face = "plain", colour = "black", size = 13),  # 设置图例文本样式
      legend.title = element_text(face = "plain", colour = "black", size = 13),  # 设置图例标题样式
      panel.grid.major = element_blank(),  # 隐藏主要网格线
      panel.grid.minor = element_blank(),  # 隐藏次要网格线
      axis.text.x = element_text(angle = 90, hjust = 0.5))  # 居中对齐
  print(p1)
  dev.off()
}

metascape2barplot_onlyUp <- function(upResFile,pdfDir,colors=c('#4497E5','#E66868')){
  library(readxl)
  library(stringr)
  library(ggplot2)
  
  ## 处理上游结果
  #upResFile <- '/Users/lin/Desktop/backup/project/shaoTong/OA/result/差异/差异_SW-M1-HTvsWT/up_ALL/metascape_result.xlsx'
  up <- read_excel(upResFile,sheet = 2)
  loc <- str_split(up$GroupID,'_',simplify = T)[,2]=='Summary'
  up <- up[loc,]
  up$color <- rep('Up',nrow(up))
  up$LogP <- -1*up$LogP
  
  ## 合并数据
  barData <- up
  
  
  ## 如果有重复
  if(sum(duplicated(barData$Description))){
    loc<-duplicated(barData$Description)
    barData$Description[loc] <- paste0(barData$Description[loc],'_2')
  }
  barData$Description <- factor(barData$Description,levels = barData$Description)
  
  
  pdf(pdfDir,width = 15)
  p1 <- ggplot(barData,aes(x = Description,y = LogP))+
    geom_col(aes(fill=color))+
    theme(axis.text.x = element_text(angle = 90))+ ## 横坐标标签调整角度
    ylab('log(p value)')+
    coord_flip() + 
    scale_fill_manual(values = colors) +
    theme(#legend.position = 'none',
      panel.background = element_rect(fill = "white", colour = "black", size = 0.25),  # 设置面板背景
      plot.title = element_text(hjust = 0.5),  # 设置标题水平居中对齐
      axis.line = element_line(colour = "black", linewidth = 0.25),  # 设置坐标轴线样式
      axis.title = element_text(size = 13, face = "plain", color = "black"),  # 设置坐标轴标题样式
      axis.text = element_text(size = 13, face = "plain", color = "black"),  # 设置坐标轴标签样式
      legend.text = element_text(face = "plain", colour = "black", size = 13),  # 设置图例文本样式
      legend.title = element_text(face = "plain", colour = "black", size = 13),  # 设置图例标题样式
      panel.grid.major = element_blank(),  # 隐藏主要网格线
      panel.grid.minor = element_blank(),  # 隐藏次要网格线
      axis.text.x = element_text(angle = 0, hjust = 0.5))  # 居中对齐
  print(p1)
  dev.off()
}


readFiles <- function(countFileDir,barcodeFileDir){
  
  ## countFileDir，barcodeFileDir需要读取count和
  ## barcode文件的地址位置
  
  count <- as.data.frame(fread(countFileDir))
  genesName <- count[,1]
  count <- count[,-1]
  row.names(count) <- genesName
  
  barcode <- as.data.frame(fread(barcodeFileDir,header = F))
  colnames(barcode) <- 'barcode'
  
  loc <- match(barcode$barcode,colnames(count))
  naLoc <- which(is.na(loc))
  if (length(naLoc)){
    print('样本中似乎没有这些barcode：')
    print(barcode$barcode[naLoc])
    print(naLoc)
  }
  
  return(list(count,barcode))
  
}

sortByBarcodeFile <- function(Count,Barcode,sortVector,n=8){
  ## Count,定量变量
  ## Barcode Barcode变量
  ## sortVector 根据barcode文件顺序对应的实验组名称向量
  ## n 每个时间节点的重复样本数量
  
  ## 根据barcode对样本进行分组，默认每组8个样本
  ## 对counts定量数据根据分组命名
  ## 后缀为组名_1到组名_8
  ## 按照barcode文件中barcode顺序对应的时间节点进行命名
  sampleName <- c()
  for (i in 1:length(sortVector)){
    sampleName <- c(sampleName,paste0(sortVector[i],'_',1:n))
  }
  
  Barcode$sampleName <- sampleName
  
  
  ## 根据barcode匹配对应的sampleName
  loc <- match(Barcode$barcode,colnames(Count))
  
  ## 如果有在定量文件中，没有被定量的barcode
  ## 则去掉这个barcode
  if (sum(is.na(loc))){
    loc <- na.omit(loc)
    ## 对定量文件重新排序和命名
    Count <- Count[,loc]
    loc2 <- match(colnames(Count),Barcode$barcode)
    colnames(Count) <- Barcode$sampleName[loc2]
  }
  else{
    ## 对定量文件重新排序和命名
    Count <- Count[,loc]
    colnames(Count) <- Barcode$sampleName
  }
  
  
  return(Count)
  
}


## 定义对细胞亚群进行分群的函数
## 输入一个Seurat对象
## 输入一个resolution向量
## 输入一个结果保存的文件夹
process_seurat_subset <- function(seuratOb,seurat_subset, result_dir,SeuratObRDS_dir, resolution_vector, npcs = 50, nfeatures = nfeatures, autoSelectDims = TRUE, dimShold = 0.1, selectPcaDim = NULL,pdf_height = 10, pdf_width = 13) {
  
  ## seuratOb，Suerat对象
  ## seurat_subset，Suerat对象
  ## result_dir，结果保存的文件夹
  ## SeuratObRDS_dir，保存RDS的文件夹
  ## resolution_vector，分辨率向量
  ## npcs，pca的维度
  ## nfeatures，寻找高变异基因的数量
  ## autoSelectDims，是否自动选择pca维度
  ## dimShold，自动选择pca维度的阈值
  ## pdf_height，保存的图片的高度
  ## pdf_width，保存的图片的宽度
  ## selectPcaDim 手动选择的pca维度

  ## 保存RDS的文件夹
  


  ## 创建结果目录
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  ## 标准化
  #seurat_subset <- NormalizeData(seurat_subset)
  ## 寻找高变异基因，并绘制火山图
  seurat_subset <- FindVariableFeatures(seurat_subset, selection.method = "vst", nfeatures = nfeatures)
  top10 <- head(VariableFeatures(seurat_subset), 10)
  pdf(paste0(result_dir, '高变异基因火山图_Fib.pdf'), width = 15, height = pdf_height)  
  plot1 <- VariableFeaturePlot(seurat_subset)    
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # 调整图片中所有字体的大小
  # 增加横坐标标签到图片的距离
  print(plot1 + plot2 & themeSet & theme(axis.text.x = element_text(margin = margin(t = 15))))
  dev.off()
  
  ## 归一化
  all.genes <- rownames(seurat_subset)
  seurat_subset <- ScaleData(seurat_subset, features = VariableFeatures(object = seurat_subset))
  gc()
  
  ## PCA
  seurat_subset <- RunPCA(seurat_subset, features = VariableFeatures(object = seurat_subset), npcs = npcs)
  print(seurat_subset[["pca"]], dims = 1:5, nfeatures = 20)
  
  ## 删除seurat_subset对象的scale.data 数据
  #seurat_subset <- DietSeurat(seurat_subset, assays = "RNA", dimreducs = "pca", scale.data = F)
  
  ## 绘制endo图
  pdf(paste0(result_dir, 'endogram_Fib.pdf'), width = pdf_width, height = pdf_height)
  print(ElbowPlot(seurat_subset, ndims = 50) & themeSet)
  dev.off()

  ## 自动选择pca维度
  if (autoSelectDims) {
    n <- length(seurat_subset@reductions[["pca"]]@stdev)
    a <- seurat_subset@reductions[["pca"]]@stdev[1:(n - 1)] - seurat_subset@reductions[["pca"]]@stdev[2:n]
    selectPcaDim <- which(a < dimShold)[1]
  }
  selectPcaDim_raw <- selectPcaDim
  print(paste0('选择的PCA维度是：', selectPcaDim))
  
  seurat_subset$group
  ## for循环resolution调用umapAndZhuShifun进行分群
  seuratList_subset <- list()
  for (resolution in resolution_vector) {
    
    ## 调用umapAndZhuShifun进行分群,并返回分群后的对象
    seuratDev <- umapAndZhuShifun(SeuratOb = seurat_subset,
                     selectPcaDim = selectPcaDim,
                     resolution = resolution,
                     resultDir = result_dir,
                     ifZhuShi = F,
                     ifplotsample = T,
                     pdf_height = pdf_height / 2,
                     pdf_width = pdf_width / 2) 


      ## 计算分群占全部细胞的占比
      cell_data <- seuratDev@meta.data %>%
      group_by(group, seurat_clusters) %>%
      summarise(count = n()) %>%
      mutate(proportion = count / ncol(SeuratOb)) %>% # 计算比例
      ungroup()

    ## 绘制分群占全部细胞的占比图
    pdf(paste0(result_dir,'分群结果_dim-',selectPcaDim,'_resolution-',resolution,'/分群占全部细胞的堆砌条形图.pdf'), width = 8, height = 8)
    p1 <- ggplot(cell_data, aes(x = group, y = proportion, fill = factor(seurat_clusters))) +
      geom_bar(stat = "identity") +
      theme(text = element_text(size=axis_text_size), legend.text = element_text(size=legend_text_size)) +
      labs(x = "Group", y = "Proportion", fill = "Cluster") +
      scale_fill_manual(values = col)+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1)) & # 使x轴标签倾斜，便于显示
      themeSet

    print(p1)
    dev.off()

    seuratList_subset[[as.character(resolution)]] <- seuratDev

  }


  
  return(seuratList_subset)
}


## 自定义cellchat流程的函数
cellchat_analysis <- function(data, pdfDir,rdsDir,species='mouse',DB='ALL',p = 10,minCell=10){
  
  ## 创建文件夹
  dir.create(pdfDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rdsDir, recursive = TRUE, showWarnings = FALSE)

  ## 选择人类或者小鼠数据库
  if(species=='human'){
   CellChatDB <- CellChatDB.human
  }
  else if(species=='mouse'){
    CellChatDB <- CellChatDB.mouse
  }
  
  
  ## 设置cellchat的数据库
  if(DB=='ALL'){   ## 全部数据库
    CellChatDB.use <- CellChatDB
  }else if(DB=='noNonProtein'){   ## 排除非蛋白信号的数据库
    CellChatDB.use <- subsetDB(CellChatDB)
  }## 否则通过参数提供指定数据库
  else {
    CellChatDB.use <- DB
  }

  ## 展示数据库信息
  pdf(paste0(pdfDir,'CellChatDB信息.pdf'),width=pdf_width,height=pdf_height)
  showDatabaseCategory(CellChatDB) & themeSet
  dev.off()

  ## 保存数据库信息
  CellChatDB_interaction <- dplyr::glimpse(CellChatDB$interaction)
  write.csv(CellChatDB_interaction,file=paste0(pdfDir,'CellChatDB_interaction.csv'))

  ##创建cellchat对象
  cellChat <- createCellChat(object = data, group.by = "ident", assay = "RNA")
  cellChat@DB <- CellChatDB.use
  cellChat <- subsetData(cellChat)

  ## 设置平行计算参数，会建立p个进程分别进行运算
  future::plan("multisession", workers = p)

  ## 识别过表达基因
  cellChat <- identifyOverExpressedGenes(cellChat)

  ## 识别过表达的互作
  cellChat <- identifyOverExpressedInteractions(cellChat)
  
  cellChat <- computeCommunProb(cellChat, type = "triMean")
  ## 保存为RDS对象
  saveRDS(cellChat,paste0(rdsDir,'chellat_computeCommunProb_res.RDS'))
  
  
  ## 输出运行时的参数
  print(cellChat@options$parameter)
  
  ## 过滤细胞
  cellChat <- filterCommunication(cellChat, min.cells = minCell)
  ## 获取细胞通讯结果
  df.net <- subsetCommunication(cellChat)
  ## 保存结果为csv文件
  write.csv(df.net,file=paste0(pdfDir,'cellchat细胞互作结果.csv'))

  ## 计算通路间的互作概率
  cellChat <- computeCommunProbPathway(cellChat)
  ## 聚合网络
  cellChat <- aggregateNet(cellChat)
  ## 保存为RDS对象
  saveRDS(cellChat,paste0(rdsDir,'cellchat_aggregateNet_res.RDS'))
  
  groupSize <- as.numeric(table(cellChat@idents))

  pdf(paste0(pdfDir,'cellchat_网络图.pdf'),width=pdf_width*2,height=pdf_height)
  par(mfrow = c(1,2))
  netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()

  ## 绘制每种细胞类型的互作网络图
  mat <- cellChat@net$weight
  cellTypes <- row.names(mat)
  ## 将cellTypes中的/替换为-
  cellTypes <- stringr::str_replace(cellTypes,'/','_')
  for (i in 1:nrow(mat)) {
    cellType <- cellTypes[i]
    pdf(paste0(pdfDir,cellType,'_互作网络图.pdf'),width=pdf_width,height=pdf_height)
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
  }

  return(cellChat)
}


run_cellchat_1 <- function(data,pdfDir,sourceSet,species = 'mouse',ifAllDB = TRUE,search=NULL,key=NULL,mem=15,workers=30,min.cells=10){

  ## 参数说明
  ## data: Seurat对象，用于cellChat分析ann
  ## pdfDir: 输出pdf的目录
  ## sourceSet: 输出图片的配置，包括ggplot的主题和图片的大小
  ## species: 物种信息，用于选择cellChat数据库
  ## ifAllDB: 是否使用全部数据库，默认为TRUE
  ## search: 搜索关键词，用于筛选数据库
  ## key: 关键词，用于筛选数据库
  ## mem: 内存大小，默认为15GB

  ## 加载cellchat配置
  source(sourceSet)

  ## 创建cellchat对象
  cellChat <- createCellChat(object = data, group.by = "ident", assay = "RNA")

  ## 选择cellChat配受体数据库
  if(species == 'mouse'){
    CellChatDB <- CellChatDB.mouse
  }else if(species == 'human'){
    CellChatDB <- CellChatDB.human
  }

  ## 创建输出目录
  dir.create(pdfDir,recursive = TRUE)

  ## 展示数据库信息
  pdf(paste0(pdfDir,'CellChatDB信息.pdf'),width=pdf_width,height=pdf_height)
  print(showDatabaseCategory(CellChatDB) & themeSet)
  dev.off()

  ## 展示数据库结构
  CellChatDB_interaction <- dplyr::glimpse(CellChatDB$interaction)
  ## 保存为csv文件
  write.csv(CellChatDB_interaction,file=paste0(pdfDir,'CellChatDB_interaction.csv'))


  ## 选择使用的数据库
  if(ifAllDB){
    CellChatDB.use <- CellChatDB
  }else{
    ## 检查是否有搜索关键词和关键词
    if(!is.null(search) & !is.null(key)){
      CellChatDB.use <- subsetDB(CellChatDB, search = list(search), key = key)
    }
    ## 否则报错
    else{
      stop('请输入有效搜索关键词和关键词')
    }
  }
  cellChat@DB <- CellChatDB.use

                         
  cellChat <- subsetData(cellChat)   

  ## 设置平行计算参数，会建立n(默认30)个进程分别进行运算
  future::plan("multisession", workers = workers)                   

  ## 识别过表达基因
  cellChat <- identifyOverExpressedGenes(cellChat)


  ## 因为下一步对全局变量大小有要求
  ## 增加全局变量大小限制
  options(future.globals.maxSize = mem * 1024^3) 

  ## 识别过表达的互作
  cellChat <- identifyOverExpressedInteractions(cellChat)


  ptm = Sys.time()
  ## type可选：triMean和truncatedMean
  ## triMean产生更少互作，但是更稳健
  cellChat <- computeCommunProb(cellChat, type = "triMean")
  Sys.time() - ptm 

  ## 运行时间长，保存关键步骤结果
  saveRDS(cellChat,file=paste0(pdfDir,'cellChat.RDS'))
  ## 查看运行的各种参数
  print(cellChat@options$parameter)

  ## 默认过滤掉细胞数量少于10个细胞的细胞群
  cellChat <- filterCommunication(cellChat, min.cells = min.cells)

  ## 获取细胞通讯结果
  df.net <- subsetCommunication(cellChat)
  ## 保存结果为csv文件
  write.csv(df.net,file=paste0(pdfDir,'cellchat细胞互作结果.csv'))

  ## 推断细胞通讯相关的功能通路
  cellChat <- computeCommunProbPathway(cellChat)

  ## 整合细胞两两间的互作信息
  cellChat <- aggregateNet(cellChat)
  ## 保存为RDS对象
  saveRDS(cellChat,paste0(pdfDir,'cellchat_aggregateNet_res.RDS'))

  groupSize <- as.numeric(table(cellChat@idents))

  pdf(paste0(pdfDir,'cellchat_网络图.pdf'),width=pdf_width*2,height=pdf_height)
  par(mfrow = c(1,2))
  netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()

  library(stringr)
  ## 单独绘制每种细胞类型的互作网络图
  mat <- cellChat@net$weight
  cellTypes <- row.names(mat)
  ## 将cellTypes中的/替换为-
  cellTypes <- str_replace(cellTypes,'/','_')
  for (i in 1:nrow(mat)) {
    cellType <- cellTypes[i]
    pdf(paste0(pdfDir,cellType,'_互作网络图.pdf'),width=pdf_width,height=pdf_height)
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    print(netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]))
    dev.off()
  }

  ## 推测模式数量
  library(NMF)
  library(ggalluvial)

  pdf(paste0(pdfDir,'cellchat_推测模式数量_outgoing.pdf'),width=pdf_width,height=pdf_height)
  print(selectK(cellChat, pattern = "outgoing"))
  dev.off()  

  pdf(paste0(pdfDir,'cellchat_推测模式数量_incoming.pdf'),width=pdf_width,height=pdf_height)
  print(selectK(cellChat, pattern = "incoming"))
  dev.off() 
  return(cellChat)
}




run_cellchat_2 <- function(cellChat,pdfDir,nPatterns,sourceSet,pattern='outgoing'){
  ## 参数说明
  ## cellChat: cellChat对象，run_cellchat_1的返回值 
  ## pdfDir: 输出pdf的目录
  ## nPatterns: 推测的模式数量  
  ## sourceSet: 输出图片的配置，包括ggplot的主题和图片的大小
  ## pattern: 推测的方向，可选：outgoing和incoming
  
  ## 加载cellchat配置
  source(sourceSet)

  ## 推测模式数量
  library(NMF)
  library(ggalluvial)

  pdf(paste0(pdfDir,paste0('cellchat_推测模式数量_',pattern,'.pdf')),width=pdf_width,height=pdf_height)
  print(selectK(cellChat, pattern = pattern))
  dev.off()

  pdf(paste0(pdfDir,paste0('cellchat_',pattern,'_heatmap.pdf')))
  cellChat <- identifyCommunicationPatterns(cellChat, pattern = pattern, k = nPatterns)
  dev.off()

  pdf(paste0(pdfDir,paste0('cellchat_',pattern,'_reverPlot.pdf')),width=pdf_width*2,height=pdf_height)
  print(netAnalysis_river(cellChat, pattern = pattern) & themeSet)
  dev.off()
  
  ## 绘制点图
  pdf(paste0(pdfDir,paste0('cellchat_',pattern,'_dotPlot.pdf')),width=pdf_width*2.5,height=pdf_height)
  print(netAnalysis_dot(cellChat, pattern = pattern,dot.size=c(3,8)) + theme(axis.text.x = element_text(angle = 90),
    panel.grid.major=element_line(size=0.5,linetype='solid',colour="gray90"),
    panel.grid.minor=element_line(size=0.5,linetype='solid',colour="gray90")
  ))
  dev.off()
  return(cellChat)
}

## 自定义计算相似行的函数



## 自定义绘制每个分组的umap图的函数
plotGroupUmap <- function(SeuratOb,groupName,pdfDir,col,pt.size = 1,label = T,ifSample=F,nSamples=NULL){
  ## 参数解释
  ## SeuratOb: Seurat对象
  ## groupName: 分组名称
  ## pdfDir: pdf存储位置
  ## col: 颜色
  ## pt.size: 点大小
  ## label: 是否显示标签
  ## ifSample: 是否对分组进行随机下采样
  ## nSamples: 随机下采样的细胞数量
  
  ## 取出分组细胞
  SeuratOb_group <- subset(x = SeuratOb,subset = group == groupName)

  ## @ date 2025-04-16
  ## 新增功能，将所有细胞类型按照原本的比例随机采样到指定的细胞数量
  ## 目的是为了让每个分组的细胞数量一致
  ## 这样做的目的是为了展示umap图的时候，可以更直观地看到每个分组的细胞类型分布的变化

  ## 当需要对分组进行随机下采样时
  if (ifSample){
    ## 取出每个细胞类型的细胞名称
    SeuratOb_group$cellName <- colnames(SeuratOb_group)
    
    SeuratOb_group_meta <- as.data.table(SeuratOb_group@meta.data)

  
    ## 取出每个细胞类型的细胞数量
    cellType_count <- table(SeuratOb_group$cellType)
    ## 计算每个细胞类型的细胞数量占总细胞数量的比例
    cellType_proportion <- cellType_count / sum(cellType_count)
    nSamples_per_type <- round(cellType_proportion * nSamples)

    ## 随机采样每个细胞类型的细胞
    ## 注意：这里的sampled_cells是一个data.table对象，不是Seurat对象
    sampled_cells <- SeuratOb_group_meta[, .SD[sample(.N, size = min(.N, nSamples_per_type[cellType]))], by = cellType]
        
    ## 输出取样前后每个细胞类型的细胞数量
    print(cellType_count)
    print(nSamples_per_type)

    ## 取出每个细胞类型的细胞名称
    sampled_cells_name <- sampled_cells$cellName

    ## 取出采样后的Seurat对象
    SeuratOb_group <- subset(x = SeuratOb_group,cells = sampled_cells_name)

  }
  
  print(levels(SeuratOb_group$cellType))
  print(dim(SeuratOb_group))
  

  pdf(paste0(pdfDir,groupName,'_聚类结果.pdf'),width=pdf_width,height=pdf_height)
  p1 <- DimPlot(SeuratOb_group,cols = col,pt.size = pt.size,label = label) & themeSet
  print(p1)
  dev.off()
}

