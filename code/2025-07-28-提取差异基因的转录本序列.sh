## @date 2025-07-28
## @author lin
## @description 提取差异基因的转录本序列,然后接着
## 在SIREs数据库中预测这些转录本序列中是否含有IRE元件


conda activate hicseq

## 定义gtf文件的位置
gtfFile="/media/ssd/sdb1/data/ljh/ref/cellranger/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf"

## 定义参考基因组的位置信息
faFiles="/media/ssd/sdb1/data/ljh/software/ref/mm10_CHIATAC/mm10/bwa/mm10.fa"

## TODO 引入差异基因列表

## check 数据
head $gtfFile
head $faFiles

# 提取GTF中的转录本信息
grep -wFf gene_list.txt ${gtfFile} | awk '$3 == "transcript" {print $1, $4-1, $5, $10, ".", $7}' | sed 's/"//g' > transcripts.bed

# 使用bedtools提取参考基因组中的转录本序列
bedtools getfasta -fi reference_genome.fa -bed transcripts.bed -fo transcripts.fasta

# 反向互补负链的序列
awk '/^>/ {if (seq) print seq; print $0; seq=""; next} {seq = seq $0} END {print seq}' transcripts.fasta | \
while read header; do
  read seq
  if [[ $header == *"-"* ]]; then
    echo "$header"
    echo "$seq" | rev | tr 'ATGC' 'TACG'
  else
    echo "$header"
    echo "$seq"
  fi
done > transcripts_with_reverse_complement.fasta
