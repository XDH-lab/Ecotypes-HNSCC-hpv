####得分
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggsci)

#细胞stat：

cell_types <- list.files(path = "/lkn_lab/znh/12月课题实验/实验/验证集/bulk/TCGA/exp", full.names = FALSE)
cell_types<-cell_types[-which(cell_types=="Ecotypes")]
# 初始化一个空的数据框来存储结果

gene_list<-list()
# 循环遍历每个细胞类型
for (cell_type in cell_types) {
  # 读取每个细胞类型的 state_assignment.txt 文件
  state_assignment_path <- paste0("/lkn_lab/znh/12月课题实验/实验/训练集/结果7//", cell_type)
  #files <- list.files(path = state_assignment_path, full.names = TRUE)
  ge <- read.table(paste0(state_assignment_path,"/gene_info.txt"), header = TRUE)
  ge$stat<-paste(cell_type,ge$Stat,sep="_")
  gene_lists <- split(ge$Gene, ge$stat)
  gene_list<-c(gene_list,gene_lists)
}


GSM8633891_z<-readRDS("/lkn_lab/znh/12月课题实验/实验/验证集/空间转录组2/不合并/GSM8633891.rds")

##看一个样本的结果 
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/空间转录组2/不合并/相关性复现美化/91")

combined_data<-GSM8633891_z
combined_data <-AddModuleScore(combined_data,features = gene_list,name = names(gene_list))
head(combined_data@meta.data)

score<-combined_data@meta.data[,10:ncol(combined_data@meta.data)]
head(score)
colnames(score)<-names(gene_list)

class<-read.table("/lkn_lab/znh/12月课题实验/实验/训练集/结果7/Ecotypes/ecotypes.txt",header=T)
head(class)

score1<-score[,class$ID]
head(score1)
spearman_corr_matrix <- cor(score1, method = "spearman")


library(corrplot)
library(reshape2)

ha <- HeatmapAnnotation(df = class[, c("Ecotype","State","Cell.type")],
                        name = "Cell Type",
                        col = list(Ecotype = c( 'E1'="#D6372E",'E2'="#FADD4B",'E3'="#70B460",'E4'="#E690C1",
              'E5'="#985EA8",'E6'="#A3A3A3"),
			  Cell.type=c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono_Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF"),
			  State=c("S01" = "#D4256B", "S02" = "#3B95AE", "S03" = "#5CA84B", "S04" = "#754B9A", "S05" = "#FF8000")))
le <- rowAnnotation(df = class[, c("Ecotype","State","Cell.type")],
                        name = "Cell Type",
                        col = list(Ecotype = c( 'E1'="#D6372E",'E2'="#FADD4B",'E3'="#70B460",'E4'="#E690C1",
              'E5'="#985EA8",'E6'="#A3A3A3"),
			  Cell.type=c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono_Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF"),
			  State=c("S01" = "#D4256B", "S02" = "#3B95AE", "S03" = "#5CA84B", "S04" = "#754B9A", "S05" = "#FF8000")))

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/空间转录组2/不合并/相关性")
pdf("特征相关性热图.pdf",width=10)
heatmap_simple(spearman_corr_matrix, 
               top_annotation = ha, # 如果没有顶部注释，可以设置为NULL
               left_annotation = le, # 如果没有左侧注释，可以设置为NULL
               width = unit(5, "in"), 
               height = unit(3, "in"), 
               legend_name = "Correlation",
               color_range = c(seq(0, quantile(spearman_corr_matrix, .8, na.rm = T), length.out = 8)), 
			   color_palette = c("gray", brewer.pal(8, "YlGnBu")), 
			   raster_quality = 20)
dev.off()
spearman_corr_matrix1<-spearman_corr_matrix
grp<-class[,c("ID","Ecotype")]
# 加载必要的包
library(dplyr)


# 根据Ecotype分组
groups <- split(grp$ID, grp$Ecotype)

# 初始化列表来存储结果
intra_correlation_matrices <- list()
inter_correlation_matrices <- list()

# 遍历每个Ecotype组
for (ecotype in names(groups)) {
  # 提取当前Ecotype的ID
  ids <- groups[[ecotype]]
  
  # 提取相同Ecotype的相关性矩阵
  intra_matrix <- spearman_corr_matrix1[ids, ids]
  
  # 存储相同Ecotype的相关性矩阵
  intra_correlation_matrices[[ecotype]] <- intra_matrix
  
  # 计算不同Ecotype的相关性矩阵
  # 首先获取所有其他Ecotype的ID
  other_ids <- setdiff(colnames(spearman_corr_matrix1), ids)
  
  # 提取不同Ecotype的相关性矩阵
  inter_matrix <- spearman_corr_matrix1[ids, other_ids]
  
  # 存储不同Ecotype的相关性矩阵
  inter_correlation_matrices[[ecotype]] <- inter_matrix
}

# 输出相同Ecotype对应ID两两相关性结果
print("Intra-Ecotype Correlation Matrices:")
print(intra_correlation_matrices)

# 输出不同Ecotype对应ID两两相关性结果
print("Inter-Ecotype Correlation Matrices:")
print(inter_correlation_matrices)
# 定义一个函数来提取相关性值和来源
extract_correlations <- function(matrices, type) {
  # 初始化一个空列表来存储提取的数据
  extracted_data <- list()
  
  # 遍历每个矩阵
  for (name in names(matrices)) {
    matrix <- matrices[[name]]
    
    # 将矩阵转换为长格式数据框
    long_matrix <- as.data.frame(as.table(matrix))
    
    # 为数据框添加一个新列，表示来源类型
    long_matrix$Source <- type
    
    # 将数据框添加到列表中
    extracted_data[[name]] <- long_matrix
  }
  
  # 将列表中的数据框合并为一个数据框
  do.call(rbind, extracted_data)
}

# 提取相同Ecotype的相关性值
intra_data <- extract_correlations(intra_correlation_matrices, "Same Ecotype")

# 提取不同Ecotype的相关性值
inter_data <- extract_correlations(inter_correlation_matrices, "Different Ecotype")

# 合并相同和不同Ecotype的数据
combined_data <- rbind(intra_data, inter_data)

# 查看合并后的数据
print(combined_data)
combined_data$Source <- factor(combined_data$Source, levels = c("Same Ecotype", "Different Ecotype"))
my_comparisons <- list( c("Same Ecotype", "Different Ecotype"))


pdf("vlo-91不处理负相关不看p值.pdf",width=4,height=4)
  ggviolin(combined_data, x = 'Source', y = 'Freq', fill = 'Source',  
               palette = c("#EB7E60","#7AC3DF"),  
               add = 'boxplot', add.params = list(fill = "white")) +
  labs(title = "Cell state spatial correlation", x = "Same ecotype", y = "Spearman correlation") +  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 'wilcox.test')  
dev.off()
