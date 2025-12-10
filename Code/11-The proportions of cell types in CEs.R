#添加生态型-TCGA样本的多特征图
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/514添加的图片/训练集-eco/热图")
rm(list=ls())
sc<-readRDS("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/转录因子分析/eco_sc.rds")
head(sc)
hpv_df<-sc@meta.data[,c(13,4)]
# 读取临床数据

hpv_table <- table(sc@meta.data[,c(13,4)])
hpv_df <- as.data.frame(hpv_table)

# 重命名列以便于理解
colnames(hpv_df) <- c("MaxEcotype", "HPV", "Count")
total_negative <- sum(hpv_df$Count[hpv_df$HPV == "HPV-"])
total_positive <- sum(hpv_df$Count[hpv_df$HPV == "HPV+"])

hpv_df$p<-hpv_df$Count
hpv_df$p[which(hpv_df$HPV=="HPV-")]<-hpv_df$Count[which(hpv_df$HPV=="HPV-")]/total_negative
hpv_df$p[which(hpv_df$HPV=="HPV+")]<-hpv_df$Count[which(hpv_df$HPV=="HPV+")]/total_positive

# 
library(tidyr)
library(ggplot2)
ls<-hpv_df
hpv_df<-hpv_df[,c(1,2,4)]
hpv_wide <- pivot_wider(hpv_df, names_from = HPV, values_from = p)
hpv_wide
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/514添加的图片/训练集-eco/热图")
n_sorted<-as.data.frame(hpv_wide)
rownames(n_sorted)<-n_sorted$MaxEcotype
n_sorted<-n_sorted[,-1]

# 保存图表
pdf("HPV_Proportion_Heatmap-g.pdf",width=2,height=4)
print(p)
dev.off()

p <- pheatmap(n_sorted, show_colnames = T,cluster_cols=F,
         border="white",#scale = "row",
         #color = c(colorRampPalette(colors = c("#4459CB","white"))(100/2),colorRampPalette(colors = c("white","#B1001F"))(100/2)),
		 #color = colorRampPalette(colors = c("white","darkgreen"))(100/2),
         cluster_rows = FALSE)

# 保存图表
pdf("HPV_Proportion_Heatmap-y.pdf",width=2,height=4)
print(p)
dev.off()



###############
##细胞类型数量百分比
library(dplyr)

# 读取数据
df <- as.data.frame(sc@meta.data)

# 统计不同 Ecotype 对应不同 recluster 的数量
count_df <- df %>%
  group_by(Ecotype, recluster) %>%
  summarise(n = n())

# 计算每个 recluster 类型的总数
total_df <- df %>%
  group_by(recluster) %>%
  summarise(total_n = n())

# 计算百分比
percentage_df <- count_df %>%
  left_join(total_df, by = "recluster") %>%
  mutate(percentage = n / total_n)

# 输出结果
print(percentage_df)

hpv_df<-as.data.frame(percentage_df[,c(1,2,5)])
hpv_wide <- pivot_wider(hpv_df, names_from = recluster, values_from = percentage)
hpv_wide
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/514添加的图片/训练集-eco/细胞分布")
n_sorted<-as.data.frame(hpv_wide)
rownames(n_sorted)<-n_sorted$Ecotype
n_sorted<-n_sorted[,-1]
which(is.na(n_sorted))
n_sorted[is.na(n_sorted)] <- 0

#bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
p <- pheatmap(n_sorted, show_colnames = T,cluster_cols=F,
         border="white",#scale = "row",
         #color = c(colorRampPalette(colors = c("#4459CB","white"))(100/2),colorRampPalette(colors = c("white","#B1001F"))(100/2)),
		 color = colorRampPalette(colors = c("white","#B1001F"))(100/2),
         cluster_rows = FALSE)

# 显示图表

# 保存图表
pdf("HPV_Proportion_Heatmap.pdf",width=6,height=4)
print(p)
dev.off()
