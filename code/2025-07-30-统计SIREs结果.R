## @date 20250730
## @description 这个脚本是用来统计SIREs网站 输出的gff格式的结果的

# conda activate RNAseq
# R

## SIREsDir="/media/ssd/sdb1/data/ljh/TFR1/result/20250730_SIREs_result/"
## tail +6 ${SIREsDir}/HSC_UP.gff > ${SIREsDir}/HSC_UP_tail.gff
## tail +6 ${SIREsDir}/HSC_UP2.gff > ${SIREsDir}/HSC_UP2_tail.gff
## cat ${SIREsDir}/HSC_UP_tail.gff ${SIREsDir}/HSC_UP2_tail.gff > ${SIREsDir}/HSC_UP.gff

## 引用包
library(data.table)
library(rtracklayer)

## 确定目录
SIREsDir <- "/media/ssd/sdb1/data/ljh/TFR1/result/20250730_SIREs_result/"
setwd(SIREsDir)

## 定义处理SIREs结果的函数
processSIREsResult <- function(gffFile, outputDir) {
  ## 读取GFF文件
  df <- import(gffFile, format = "gff")
  df <- as.data.table(df)
  
  ## 去掉字符型列中的引号
  df[] <- lapply(df, function(x) {
    if(is.character(x)) {
      gsub('"', '', x)
    } else {
      x
    }
  })
  
  ## 输出基本统计信息
  print("基本统计信息:")
  print(head(df))
  print("GU-UG数量分布:")
  print(table(df$N_GU_UG))
  print("质量分布:")
  print(table(df$QUALITY))
  print("自由能统计:")
  print(summary(as.numeric(df$free_energy)))
  print("序列位置分布:")
  print(table(df$seq_pos))
  
  ## 保存CSV结果
  outputFile <- file.path(outputDir, basename(tools::file_path_sans_ext(gffFile))) 
  write.csv(df, file = paste0(outputFile, ".csv"), row.names = FALSE)
  
  ## 处理基因名称
  df$geneName <- gsub('-.*$', '', df$seqnames)
  
  ## 提取高质量结果
  df_high <- df[QUALITY %in% c('High','High-Medium')]
  high_quality_genes <- unique(df_high$geneName)
  
  ## 返回处理后的数据框和高质量基因列表
  return(list(
    full_results = df,
    high_quality_results = df_high,
    high_quality_genes = high_quality_genes
  ))
}

## 使用函数处理数据
HSC_up_results <- processSIREsResult(
  gffFile = paste0(SIREsDir,'/HSC_UP.gff'),
  outputDir = SIREsDir
)

## 获取处理后的数据框
HSC_up_df <- HSC_up_results$full_results
HSC_up_df_high <- HSC_up_results$high_quality_results


"Nr4a3" %in% HSC_up_df$geneName
HSC_up_df[geneName == 'Nr4a3',]


## TODO 
## 1. 看SIREs的文档，看看需不需要从结果里筛选掉一些结果出去
## 2. 提取出基因list 
## 3. 取交集