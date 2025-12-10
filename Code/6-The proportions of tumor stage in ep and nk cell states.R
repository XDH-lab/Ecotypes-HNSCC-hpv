#绘制堆叠箱线图统计肿瘤分期分布情况
assi<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/结果/TCGA-HNSCC最后/Epithelial_cells/state_assignment.txt",header=T)
head(assi)
clin<-read.delim("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/clin.txt")
head(clin)

clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="Stage IVA")]<-"Stage IV"
clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="Stage IVB")]<-"Stage IV"
clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="Stage IVC")]<-"Stage IV"
table(clin$ajcc_pathologic_tumor_stage)

clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="[Discrepancy]")]<-"NA"
clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="")]<-"NA"

  
clin<-clin[,c("ID","ajcc_pathologic_tumor_stage")]
assi<-assi[,c("ID","State")]
clin$ID <- gsub("-", ".", clin$ID)
head(clin)


##排除NA情况
merged_data <- merge(clin, assi, by = "ID", all.x = TRUE)

head(merged_data)
merged_data<-na.omit(merged_data)
merged_data<-merged_data[-which(merged_data$ajcc_pathologic_tumor_stage=="NA"),]
  group_counts <- merged_data %>%
    group_by(State, ajcc_pathologic_tumor_stage) %>%
    tally() %>%
    ungroup()
	
  group_proportions <- group_counts %>%
    group_by(State) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

# 绘图group_proportions
colors <- c("Stage I" = "#56B4E9", "Stage II" = "#E69F00", "Stage III" = "#D55E00", "Stage IV" = "#0072B2")
ggplot(group_proportions, aes(x = State, y = proportion, fill = ajcc_pathologic_tumor_stage)) +
  geom_bar(stat = "identity", position = "stack") +  # 使用 position = "stack" 来堆叠柱状图
  labs(title = "Proportion of Histological Grades by Ecotype",  # 替换为具体的标题
       x = "Max Ecotype",  # 替换为具体的 x 轴标签
       y = "Proportion",
       fill = "Histological Grade") +  # 替换为具体的图例标题
  theme_minimal() +  # 选择一个主题
  scale_fill_manual(values = colors)  +  
  theme_classic()# 自定义颜色
setwd("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/肿瘤分期比例/stat/EP")
ggsave("stage_EP_percent.pdf",width=6,height=5)



#################################
assi<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/结果/TCGA-HNSCC最后/NK_cells/state_assignment.txt",header=T)
head(assi)
clin<-read.delim("/lkn_lab/znh/12月课题实验/实验/验证集/重新处理TCGA的结果/clin.txt")
head(clin)
clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="Stage IVA")]<-"Stage IV"
clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="Stage IVB")]<-"Stage IV"
clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="Stage IVC")]<-"Stage IV"
table(clin$ajcc_pathologic_tumor_stage)

clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="[Discrepancy]")]<-"NA"
clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="")]<-"NA"

# 
clin<-clin[,c("ID","ajcc_pathologic_tumor_stage")]
assi<-assi[,c("ID","State")]
clin$ID <- gsub("-", ".", clin$ID)
head(clin)
merged_data <- merge(clin, assi, by = "ID", all.x = TRUE)

head(merged_data)
merged_data<-na.omit(merged_data)
merged_data<-merged_data[-which(merged_data$ajcc_pathologic_tumor_stage=="NA"),]
  group_counts <- merged_data %>%
    group_by(State, ajcc_pathologic_tumor_stage) %>%
    tally() %>%
    ungroup()
	
  group_proportions <- group_counts %>%
    group_by(State) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

# 绘图group_proportions
colors <- c("Stage I" = "#56B4E9", "Stage II" = "#E69F00", "Stage III" = "#D55E00", "Stage IV" = "#0072B2")
ggplot(group_proportions, aes(x = State, y = proportion, fill = ajcc_pathologic_tumor_stage)) +
  geom_bar(stat = "identity", position = "stack") +  # 使用 position = "stack" 来堆叠柱状图
  labs(title = "Proportion of Histological Grades by Ecotype",  # 替换为具体的标题
       x = "Max Ecotype",  # 替换为具体的 x 轴标签
       y = "Proportion",
       fill = "Histological Grade") +  # 替换为具体的图例标题
  theme_minimal() +  # 选择一个主题
  scale_fill_manual(values = colors)  +  
  theme_classic()# 自定义颜色
ggsave("stage_nk_percent.pdf",width=6,height=5)







