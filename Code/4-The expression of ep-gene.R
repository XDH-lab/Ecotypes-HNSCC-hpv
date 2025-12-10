rm(list=ls())
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot)
library(ggpubr);library(beeswarm);library(presto);
library(GSEABase)

library(GSVA)

a = read.delim("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/TCGA-HNSCC最后.txt")#
rownames(a)<-a[,1]
a=a[,-1]
dim(a)
a[1:5,1:5]
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/基因表达")
lab<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/结果/TCGA-HNSCC最后/Epithelial_cells/state_assignment.txt",header=T)
scor<-t(a[c("CXCL10","MAGEA4"),])
scor<-as.data.frame(scor)
scor$ID<-rownames(scor)

# 绘制箱线图
col<-c("S01" = "#D4256B", "S02" = "#3B95AE", "S03" = "#5CA84B", "S04" = "#754B9A", "S05" = "#FF8000")

my_comparisons <- list(
  c("S03","S01"), 
  c("S03", "S02"), 
  c("S03", "S04"), c("S03", "S05")
)	
p<-ggplot(lnc, aes(x = factor(State, levels = names(col)), y = CXCL10, fill = State)) +
  geom_boxplot() +
  scale_fill_manual(values = col) +
  labs(title = "CXCL10", x = "State", y = "CXCL10 expression") +
  theme_classic()
p<-p+stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p
ggsave("TCGA-CCL10.pdf",width=3.5,height=3.5)

my_comparisons <- list(
  c("S05", "S04"),c("S05", "S03"), 
  c("S05", "S02"),c("S05","S01") 
)	
p<-ggplot(lnc, aes(x = factor(State, levels = names(col)), y = MAGEA4, fill = State)) +
  geom_boxplot() +
  scale_fill_manual(values = col) +
  labs(title = "MAGEA4 ", x = "State", y = "MAGEA4  expression") +
  theme_classic()
p<-p+stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p
ggsave("TCGA-MAGEA4.pdf",width=3.5,height=3.5)



##单细胞
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/药物敏感性/生存/基因-STAT/单细胞/ep")

sc<-readRDS("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/scdata/Epithelial_cells.rds")
data <- as.data.frame(sc@meta.data)
genes <- c("CXCL10","MAGEA4")
for (gene in genes) {
  expr_data <- FetchData(sc, vars = gene)
  data[[gene]] <- expr_data
}
my_comparisons <- list(c("S05", "S04"), c("S05", "S03"), c("S05", "S02"),c("S05", "S01"))
my_comparisons <- list(
  c("S01","S03"), 
  c("S01", "S02"), 
  c("S01", "S04"), c("S01", "S05")
)
data$stat <- factor(data$stat)
ls<-data
for (gene in genes) {
  data[[gene]] <- as.numeric(data[[gene]][,1])
  #data <- na.omit(data)  # 删除包含缺失值的行
}

my_comparisons <- list(
  c("S05", "S04"),c("S05", "S03"), 
  c("S05", "S02"),c("S05","S01") 
)	
col<-c("S01" = "#D4256B", "S02" = "#3B95AE", "S03" = "#5CA84B", "S04" = "#754B9A", "S05" = "#FF8000")

p<-ggplot(data, aes(x = factor(stat, levels = names(col)), y = MAGEA4, fill = stat)) +
  geom_boxplot() +
  scale_fill_manual(values = col) +
  labs(title = "MAGEA4 ", x = "stat", y = "MAGEA4  expression") +
  theme_classic()
p<-p+stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p
ggsave("MAGEA4-sc.pdf",width=3.5,height=3.5)

my_comparisons <- list(
  c("S05", "S04"),c("S05", "S03"), 
  c("S05", "S02"),c("S05","S01") 
)	
p<-ggplot(data, aes(x = factor(stat, levels = names(col)), y = CXCL10, fill = stat)) +
  geom_boxplot() +
  scale_fill_manual(values = col) +
  labs(title = "CXCL10 ", x = "stat", y = "CXCL10  expression") +
  theme_classic()
p<-p+stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p
ggsave("TCGA-CXCL10.pdf",width=3.5,height=3.5)
