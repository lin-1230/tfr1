## 设置随机数种子，保证结果可复现
set.seed(666)

## 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)


## 变量
nfeatures <- 2000 # 高变基因筛选数量(建议值2000-3000)
## NOTE: 该参数直接影响后续PCA和聚类效果，值过小会丢失生物学差异信息
npcs=50 # RunPCA函数的参数，用来控制最后需要生成几个pca，默认是50，在细胞数少于50（Drug-seq）的时候需要更改
autoSelectDims <- T #是否使用for循环自动挑选可能是局部最优解的dim
dimShold <- 0.1 #自动选择dim维度时的阈值

## ggplot2 图片相关设置
pt.size <- 2
axis_text_size <- 20  # 全局坐标轴字体大小
legend_text_size <- 18  # 全局图例字体大小
label_text_size <- 20  # 全局标签字体大小
title_text_size <- 20  # 全局标题字体大小

pdf_width <- 13 # PDF图形默认宽度
pdf_height <- 10 # PDF图形默认高度


## 全局颜色定义
col = c('#5668E4','#FFA9AF','#E66868','#969696','#4497E5',
        '#AFE4EA','#F8CD42','#C486CF','#CC19EF','#429342',
        '#5DFCD5','#FFECB3FF','#FCF78A','#F59D30',
        '#B59C93','#C6312C','#8994AC','#254365','#BAC0F2',
        '#E7E7E7','#A2CCEA','#835C27','#433825','#CFA670')

## 全局ggplot2主题设置
themeSet <- theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),  # 设置面板背景
        plot.title = element_text(hjust = 0.5, size = title_text_size),  # 设置标题水平居中对齐并增加字体大小
        strip.text = element_text(size = title_text_size,  # 继承标题字号
        face = "bold",          # 加粗显示
        ),  # 设置strip文本样式
        axis.line = element_line(colour = "black", linewidth = 0.25),  # 设置坐标轴线样式
        axis.title = element_text(size = (axis_text_size + 2), face = "plain", color = "black"),  # 设置坐标轴标题样式
        axis.text = element_text(size = axis_text_size, face = "plain", color = "black"),  # 设置坐标轴标签样式
        legend.text = element_text(face = "plain", colour = "black", size = (legend_text_size)),  # 设置图例文本样式
        legend.title = element_text(face = "plain", colour = "black", size = (legend_text_size+2)),  # 设置图例标题样式
        panel.grid.major = element_blank(),  # 隐藏主要网格线
        panel.grid.minor = element_blank(),  # 隐藏次要网格线
        axis.text.x = element_text(angle = 0, hjust = 0.5))  # 居中对齐

