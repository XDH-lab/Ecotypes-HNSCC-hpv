rm(list=ls())
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(reshape)

setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7/Ecotypes")

state_ecotype<-read.table("ecotypes.txt",stringsAsFactors=F,header=T)
state_ecotype<-state_ecotype[,c(1,6)]
state_ecotype$link<-1

state_ecotype<-reshape::melt(state_ecotype,id='link')
variable<-summary(state_ecotype$variable)
state_ecotype$flow<-rep(1:variable[1],length(variable))
head(state_ecotype)

####sankey
library(ggalluvial)
library(RColorBrewer)
library(randomcoloR)

col_sankey<-c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 
              'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono_Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF",
              'E1'="#D6372E",'E2'="#FADD4B",'E3'="#70B460",'E4'="#E690C1",
              'E5'="#985EA8",'E6'="#A3A3A3")#,'E7'="#B7D3E5"
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/文章细胞组成复现")

pdf("state_ecotype_sankey2.pdf",width = 6,height = 8)
ggplot(state_ecotype, aes(x = variable, y = link,
                          stratum = value, alluvium = flow, fill = value)) +
  geom_stratum() + 
  geom_flow(aes.flow = 'forward') + 
  scale_fill_manual(values = col_sankey[rank(1:44)]) + 
  guides(fill=FALSE)+
  geom_text(stat = 'stratum', infer.label = TRUE, size = 2.5) + 
  labs(x = '', y = '') + 
  theme(legend.position = 'none', panel.background = element_blank(),
        line = element_blank(), axis.text.y = element_blank())+
  scale_x_discrete(limits = c('CellType', 'Ecotype'))

dev.off()  


path <- "/lkn_lab/znh/12月课题实验/实验/训练集/结果7/"

# 获取路径下的所有文件名称
files <- list.files(path = path, full.names = FALSE, recursive = FALSE)

# 打印结果
print(files)
# 设置基础路径
base_path <- "/lkn_lab/znh/12月课题实验/实验/训练集/结果7/"

# 获取所有细胞类型文件夹
cell_types <- files[1:13]
cell_types<-cell_types[-5]
all_data <- data.frame()

# 循环读取每个文件夹下的 state_assignment.txt 文件
# 初始化一个空的数据框用于合并所有文件
all_data <- data.frame()

# 循环读取每个文件夹下的 state_assignment.txt 文件
for (cell_type in cell_types) {
  # 构建每个细胞类型文件夹的完整路径
  cell_path <- paste0(base_path, cell_type, sep = "/")

  # 检查 state_assignment.txt 文件是否存在
  if ("state_assignment.txt" %in% list.files(path = cell_path)) {
    # 读取 state_assignment.txt 文件
    state_assignment_file <- paste0(cell_path, "state_assignment.txt", sep = "")
    
    # 读取文件内容
    state_data <- read.table(state_assignment_file, header = TRUE, sep = "\t")
    
    # 添加 cell_types 列
    state_data$cell_type <- cell_type
    
    # 将数据添加到 all_data 数据框
    all_data <- rbind(all_data, state_data)
  } else {
    print(paste(cell_type, "does not contain state_assignment.txt file."))
  }
}

setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7/Ecotypes")
state_ecotype<-read.table("ecotypes.txt",stringsAsFactors=F,header=T)

ls<-all_data
all_data<-ls
head(all_data)

all_data$ID<-paste(all_data$cell_type,all_data$State,sep="_")
all_data_with_ecotype <- left_join(all_data, state_ecotype[, c("ID", "Ecotype")], by = "ID")

library(dplyr)
library(tidyr)

# all_data_with_ecotype 是已经合并好的数据框
# 计算每个 Ecotype 中每种 cell_type 的数量
count_data <- all_data_with_ecotype %>%
  group_by(Ecotype, cell_type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个 Ecotype 的总数
total_count <- count_data %>%
  group_by(Ecotype) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 合并总数数据，计算比例
proportion_data <- count_data %>%
  left_join(total_count, by = "Ecotype") %>%
  mutate(percent = count / total)
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot)
library(ggpubr);library(beeswarm);library(presto);
library(GSEABase)
mytheme <- theme_bw() + 
  theme(plot.title=element_text(size=rel(2),hjust=0.5),
        axis.title=element_text(size=rel(1.2)),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = rel(1.5), color = 'black'),
        axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.5),color = 'black'),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        panel.border=element_rect(color="black",size=1),
        axis.line=element_line(color="black",size=0.5))
col_sankey<-c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 
              'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono_Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF",
              'E1'="#D6372E",'E2'="#FADD4B",'E3'="#70B460",'E4'="#E690C1",
              'E5'="#985EA8",'E6'="#A3A3A3")#,'E7'="#B7D3E5
p = ggplot(proportion_data, aes( x = Ecotype, weight = percent, fill = cell_type))+
  geom_bar( position = "stack")+scale_fill_manual( values = col_sankey)+mytheme+labs(y = "Percent")
p
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/文章细胞组成复现")
ggsave("eco_cell_percent.pdf")