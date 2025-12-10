library(Seurat)
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot)
library(GenomicRanges);library(presto)
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)
#3种方法删除双细胞需求的R包
library(Seurat)
library(scDblFinder)
library(scds)
library(scran)

library(ggplot2)
library(dplyr)
library(tidyr)

library(pheatmap)

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/药物敏感性/生存/美化热图")
result<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/药物敏感性/生存/drg-survival.txt", header = TRUE)
# 查看结果
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/药物敏感性/生存")
sur<-read.table("drg-survival.txt",heade=T)
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/药物敏感性/生存/美化热图")
sorted_number_matrix<-sur[,c("Cisplatin_1005","Paclitaxel_1080","X5.Fluorouracil_1073","Docetaxel_1819",
                  "Afatinib_1032","Palbociclib_1054","Ribociclib_1632","Erlotinib_1168","Gefitinib_1010","Pyridostatin_2044",
				  "JAK_8517_1739",
                  "Venetoclax_1909")]
head(sorted_number_matrix)
rownames(sorted_number_matrix)<-sur$CellType_states

library(reshape2)
sorted_number_matrix_df <- as.data.frame(sorted_number_matrix, row.names = rownames(sorted_number_matrix))
head(sorted_number_matrix_df)
sorted_number_matrix_df$CellType <- rownames(sorted_number_matrix)

sorted_number_matrix_long <- reshape2::melt(sorted_number_matrix_df, 
                                  id.vars = "CellType", 
                                  variable.name = "Drug", 
                                  value.name = "Value")
dim(sorted_number_matrix_long)
head(sorted_number_matrix_long)
data<-sorted_number_matrix_long
data$Shape <- ifelse(data$Value > 0, "circle", "square")
data$Color <- with(data, ifelse(Value > 0.1, "#EFABD2",
                   ifelse(Value > 0.05, "#EA81BD",
                   ifelse(Value > 0, "#E23770",
                   ifelse(Value > -0.05, "#0D9DBF",
                   ifelse(Value > -0.1, "#2BD2D6", "#B6EFF2"))))))
head(data)	
data$Value2 <- ifelse(data$Value == 0, NA, -log10(abs(data$Value)))

ggplot(data, aes(x = CellType, y = Drug, size = Value2, shape = Shape, color = Color, fill = Color)) +
  geom_point() +  # 使用点来表示数据
  scale_size(range = c(1, 10), name = "-log10(p-value)") +  # 设置点的大小范围和名称
  scale_shape_manual(values = c(21, 22)) +  # 设置形状，21为圆形，22为方形
  scale_color_manual(values = c("#0D9DBF", "#2BD2D6", "#B6EFF2", "#E23770", "#EA81BD", "#EFABD2")) +  # 设置颜色
  scale_fill_manual(values = c("#0D9DBF", "#2BD2D6", "#B6EFF2", "#E23770", "#EA81BD", "#EFABD2"))+
  theme_classic() +  # 使用经典主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 调整x轴标签角度
        axis.text.y = element_text(angle = 0)) +  # 保持y轴标签水平
  labs(x = "Cell Type", y = "Drug")  # 设置坐标轴标签
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/药物敏感性/生存/美化热图")

# 显示图形
ggsave("heatmap_plot.pdf", width = 15, height = 5)



