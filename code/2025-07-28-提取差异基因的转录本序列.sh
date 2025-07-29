## @date 2025-07-28
## @author lin
## @description 提取差异基因的转录本序列,然后接着
## 在SIREs数据库中预测这些转录本序列中是否含有IRE元件


conda activate hicseq
cd /media/ssd/sdb1/data/ljh/TFR1/

## 定义gtf文件的位置
gtfFile="/media/ssd/sdb1/data/ljh/ref/cellranger/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf"

## 定义参考基因组的位置信息
faFiles="/media/ssd/sdb1/data/ljh/software/ref/mm10_CHIATAC/mm10/bwa/mm10.fa"

## check 数据
head -10 $gtfFile
head $faFiles


DEG_dir="/media/ssd/sdb1/data/ljh/TFR1/result/20250726/"

## 定义差异基因文件夹
DEGListDir_0p5=${DEG_dir}DEGeneList_0p5/

## ALL
ALL_upGeneList_0p5=${DEGListDir_0p5}ALL_DEG_up.csv
ALL_downGeneList_0p5=${DEGListDir_0p5}ALL_DEG_down.csv
## HSC
HSC_upGeneList_0p5=${DEGListDir_0p5}HSC_DEG_up.csv
HSC_downGeneList_0p5=${DEGListDir_0p5}HSC_DEG_down.csv
## MPP1
MPP1_upGeneList_0p5=${DEGListDir_0p5}MPP1_DEG_up.csv
MPP1_downGeneList_0p5=${DEGListDir_0p5}MPP1_DEG_down.csv
## MPP2_3
MPP2_3_upGeneList_0p5=${DEGListDir_0p5}MPP2_3_DEG_up.csv
MPP2_3_downGeneList_0p5=${DEGListDir_0p5}MPP2_3_DEG_down.csv
## cycling_MPP2_3
cMPP2_3_upGeneList_0p5=${DEGListDir_0p5}cycling_MPP2_3_DEG_up.csv
cMPP2_3_downGeneList_0p5=${DEGListDir_0p5}cycling_MPP2_3_DEG_down.csv
## MPP4
MPP4_upGeneList_0p5=${DEGListDir_0p5}MPP4_DEG_up.csv
MPP4_downGeneList_0p5=${DEGListDir_0p5}MPP4_DEG_down.csv


## 定义结果文件夹
resultDir='/media/ssd/sdb1/data/ljh/TFR1/result/20250728/'

## 从上调基因列表中提取第一列并去掉表头
cut -f 9 ${gtfFile} | head -10

## 提取差异基因的转录本位置信息

# 去掉表头，提取第一列
# 定义函数用于从GTF文件中提取转录本信息
extract_transcript_info() {
    local gene_list=$1  # 基因列表文件
    local gtf_file=$2   # GTF文件
    local outBedFiles=$3 # 输出的bed文件

    tail -n +2 ${gene_list} | cut -f1 -d "," \
    | grep -wFf - ${gtf_file} \
    | awk 'BEGIN{OFS="\t"} $3=="transcript" {
        # 从属性字段中提取transcript_name
        if (match($0, /transcript_name "[^"]+"/)) {
            transcript_name=substr($0,RSTART+17,RLENGTH-18)
            print $1, $4-1, $5, transcript_name, ".", $7
          }
    }' > ${outBedFiles}
}

mkdir -p ${resultDir}/0p5/

# 调用函数
# ALL
extract_transcript_info ${ALL_upGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_ALL_up.bed
extract_transcript_info ${ALL_downGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_ALL_down.bed
head ${resultDir}/0p5/transcripts_ALL_up.bed
head ${resultDir}/0p5/transcripts_ALL_down.bed


# HSC
extract_transcript_info ${HSC_upGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_HSC_up.bed
extract_transcript_info ${HSC_downGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_HSC_down.bed
head ${resultDir}/0p5/transcripts_HSC_up.bed
head ${resultDir}/0p5/transcripts_HSC_down.bed

# MPP1
extract_transcript_info ${MPP1_upGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_MPP1_up.bed
extract_transcript_info ${MPP1_downGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_MPP1_down.bed
head ${resultDir}/0p5/transcripts_MPP1_up.bed
head ${resultDir}/0p5/transcripts_MPP1_down.bed

# MPP2_3
extract_transcript_info ${MPP2_3_upGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_MPP2_3_up.bed
extract_transcript_info ${MPP2_3_downGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_MPP2_3_down.bed
head ${resultDir}/0p5/transcripts_MPP2_3_up.bed
head ${resultDir}/0p5/transcripts_MPP2_3_down.bed

# cycling_MPP2_3
extract_transcript_info ${cMPP2_3_upGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_cMPP2_3_up.bed
extract_transcript_info ${cMPP2_3_downGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_cMPP2_3_down.bed
head ${resultDir}/0p5/transcripts_cMPP2_3_up.bed
head ${resultDir}/0p5/transcripts_cMPP2_3_down.bed


# MPP4
extract_transcript_info ${MPP4_upGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_MPP4_up.bed
extract_transcript_info ${MPP4_downGeneList_0p5} ${gtfFile} ${resultDir}/0p5/transcripts_MPP4_down.bed
head ${resultDir}/0p5/transcripts_MPP4_up.bed
head ${resultDir}/0p5/transcripts_MPP4_down.bed



# 提取GTF中的转录本信息
# 使用bedtools提取参考基因组中的转录本序列
## ALL
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_ALL_up.bed -fo ${resultDir}/0p5/transcripts_ALL_up.fasta
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_ALL_down.bed -fo ${resultDir}/0p5/transcripts_ALL_down.fasta

## HSC
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_HSC_up.bed -fo ${resultDir}/0p5/transcripts_HSC_up.fasta
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_HSC_down.bed -fo ${resultDir}/0p5/transcripts_HSC_down.fasta

## MPP1
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_MPP1_up.bed -fo ${resultDir}/0p5/transcripts_MPP1_up.fasta
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_MPP1_down.bed -fo ${resultDir}/0p5/transcripts_MPP1_down.fasta

## MPP2_3
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_MPP2_3_up.bed -fo ${resultDir}/0p5/transcripts_MPP2_3_up.fasta
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_MPP2_3_down.bed -fo ${resultDir}/0p5/transcripts_MPP2_3_down.fasta

## cycling_MPP2_3
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_cMPP2_3_up.bed -fo ${resultDir}/0p5/transcripts_cMPP2_3_up.fasta
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_cMPP2_3_down.bed -fo ${resultDir}/0p5/transcripts_cMPP2_3_down.fasta

## MPP4
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_MPP4_up.bed -fo ${resultDir}/0p5/transcripts_MPP4_up.fasta
bedtools getfasta -name -fi ${faFiles} -bed ${resultDir}/0p5/transcripts_MPP4_down.bed -fo ${resultDir}/0p5/transcripts_MPP4_down.fasta

head -2 ${resultDir}/0p5/transcripts_MPP4_up.fasta

## 检查所有生成的fasta文件
ls -l ${resultDir}/0p5/*.fasta

wc -l ${resultDir}/0p5/*.fasta

