
#####################细胞丰度预测效果
bulk<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/bulk合并4/预后/丰度预后/生存results.txt")
head(bulk)
dim(bulk)
TCGA<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/丰度预后/生存results.txt")
head(TCGA)
dim(TCGA)

data<-TCGA[,c(1,2,4)]
head(data)
colnames(data)[3]<-"TCGA"
data$BULK<-bulk$p.value
head(data)
colors <- c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 
              'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono/Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF")
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/丰度预后一致性")
data$TCGA_log10 <- -log10(data$TCGA)
data$BULK_log10 <- -log10(data$BULK)

data$TCGA_log10 <- ifelse(TCGA$coef > 0, data$TCGA_log10, -data$TCGA_log10)
data$BULK_log10 <- ifelse(bulk$coef > 0, data$BULK_log10, -data$BULK_log10)

# 绘制散点图
# 计算相关性系数
#Spearman
correlation <- cor(data$TCGA_log10, data$BULK_log10, method = "spearman")
cor_text <- paste("Spearman Correlation: ", round(correlation, 2))

# 拟合线性模型
model <- lm(BULK_log10 ~ TCGA_log10, data = data)
model_summary <- summary(model)


coef_text <- paste("Slope: ", round(coef(model)[2], 2))
p_value_text <- paste("P-value: ", model_summary$coefficients[2, 4])

# 绘制散点图
p <- ggplot(data, aes(x = TCGA_log10, y = BULK_log10, color = CellType)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = colors) +
  labs(title = "Correlation Scatter Plot",
       x = "-log10(TCGA)",
       y = "-log10(BULK)",
       color = "Cell Type",
       caption = paste(cor_text, coef_text, p_value_text, sep = "\n")) +
  theme_classic() +
  geom_smooth(method = "lm", se = TRUE, color = "black")   # 添加线性回归线和置信区间

# 显示图表
print(p)

# 保存图表
ggsave("丰度一致性_Spearman-置信区间.pdf", plot = p, width = 5, height = 4.5)

