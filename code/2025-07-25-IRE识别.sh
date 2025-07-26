

## 允许一个错配


conda activate ATACseq
cd /media/ssd/sdb1/data/ljh/TFR1/

resDir="result/2025-07-25-IRE"
mkdir -p $resDir

## 生成IRE序列，不允许错配
seq2profile.pl CAGTGN 0 IRE_CAGTGN_0_motif > $resDir/IRE_CAGTGN_0.motif
## 生成IRE序列，允许一个错配
seq2profile.pl CAGTGN 1 IRE_CAGTGN_1_motif > $resDir/IRE_CAGTGN_1.motif
## 生成IRE序列，允许两个错配
seq2profile.pl CAGTGN 2 IRE_CAGTGN_2_motif > $resDir/IRE_CAGTGN_2.motif


head -10 $resDir/IRE_CAGTGN_0.motif
head -10 $resDir/IRE_CAGTGN_1.motif
head -10 $resDir/IRE_CAGTGN_2.motif


## 扫描全基因组，看哪些基因包含这些序列
nohup scanMotifGenomeWide.pl $resDir/IRE_CAGTGN_0.motif mm10 -bed -keepAll > $resDir/IRE_CAGTGN_0.bed 2> $resDir/IRE_CAGTGN_0.log &
nohup scanMotifGenomeWide.pl $resDir/IRE_CAGTGN_1.motif mm10 -bed -keepAll > $resDir/IRE_CAGTGN_1.bed 2> $resDir/IRE_CAGTGN_1.log &
nohup scanMotifGenomeWide.pl $resDir/IRE_CAGTGN_2.motif mm10 -bed -keepAll > $resDir/IRE_CAGTGN_2.bed 2> $resDir/IRE_CAGTGN_2.log &

tail -f $resDir/IRE_CAGTGN_0.bed

## 使用基因区域查找呢？

mm10GeneBodyBed='/media/ssd/sdc1/data/ljh/CHIATAC/data/ref/mm10_geneBody.bed'
head -2 $mm10GeneBodyBed

mkdir -p $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/
mkdir -p $resDir/homer_geneBody/res_IRE_CAGTGN_1.motif/
mkdir -p $resDir/homer_geneBody/res_IRE_CAGTGN_2.motif/

## 0个错配
nohup findMotifsGenome.pl $mm10GeneBodyBed mm10 $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/ \
                    -size given \
                    -find $resDir/IRE_CAGTGN_0.motif \
                    > $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/output.txt 2> $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/error.log &

## 1个错配
nohup findMotifsGenome.pl $mm10GeneBodyBed mm10 $resDir/homer_geneBody/res_IRE_CAGTGN_1.motif/ \
                    -size given \
                    -find $resDir/IRE_CAGTGN_1.motif \
                    > $resDir/homer_geneBody/res_IRE_CAGTGN_1.motif/output.txt 2> $resDir/homer_geneBody/res_IRE_CAGTGN_1.motif/error.log &

## 2个错配
nohup findMotifsGenome.pl $mm10GeneBodyBed mm10 $resDir/homer_geneBody/res_IRE_CAGTGN_2.motif/ \
                    -size given \
                    -find $resDir/IRE_CAGTGN_2.motif \
                    > $resDir/homer_geneBody/res_IRE_CAGTGN_2.motif/output.txt 2> $resDir/homer_geneBody/res_IRE_CAGTGN_2.motif/error.log &


## TODO 去掉GM基因，然后去重，得到基因list


## 0个错配的结果
## 从output.txt提取第一列，去掉Gm开头的基因并去重
cat $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/output.txt | \
    cut -f 1 | \
    grep -v '^Gm' | \
    sort -u > $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/IRE_CAGTGN_0_filtered_genes.txt

head $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/IRE_CAGTGN_0_filtered_genes.txt
tail $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/IRE_CAGTGN_0_filtered_genes.txt
wc -l $resDir/homer_geneBody/res_IRE_CAGTGN_0.motif/IRE_CAGTGN_0_filtered_genes.txt

## 1个错配的结果
## 从output.txt提取第一列，去掉Gm开头的基因并去重
cat $resDir/homer_geneBody/res_IRE_CAGTGN_1.motif/output.txt | \
    cut -f 1 | \
    grep -v '^Gm' | \
    sort -u > $resDir/homer_geneBody/res_IRE_CAGTGN_1.motif/IRE_CAGTGN_1_filtered_genes.txt

head $resDir/homer_geneBody/res_IRE_CAGTGN_1.motif/IRE_CAGTGN_1_filtered_genes.txt
wc -l $resDir/homer_geneBody/res_IRE_CAGTGN_1.motif/IRE_CAGTGN_1_filtered_genes.txt

## 2个错配的结果
## 从output.txt提取第一列，去掉Gm开头的基因并去重
cat $resDir/homer_geneBody/res_IRE_CAGTGN_2.motif/output.txt | \
    cut -f 1 | \
    grep -v '^Gm' | \
    sort -u > $resDir/homer_geneBody/res_IRE_CAGTGN_2.motif/IRE_CAGTGN_2_filtered_genes.txt

head $resDir/homer_geneBody/res_IRE_CAGTGN_2.motif/IRE_CAGTGN_2_filtered_genes.txt
wc -l $resDir/homer_geneBody/res_IRE_CAGTGN_2.motif/IRE_CAGTGN_2_filtered_genes.txt
