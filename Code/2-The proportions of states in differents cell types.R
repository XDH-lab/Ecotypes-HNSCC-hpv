setwd("/groups/g900008/home/zhangnihui/新课题/实验/训练GSE234933+GSE182227/修改/")
scobj<-readRDS("recluster.rds")
color<-c("S01" = "#D4256B", "S02" = "#3B95AE", "S03" = "#5CA84B", "S04" = "#754B9A", "S05" = "#FF8000")

# 设置路径
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
data<-data.frame()
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
  coun<-as.data.frame(as.matrix(table(state_data$State),,2))
  for(i in 1:nrow(coun)){
  coun$percent[i]<-coun[i,1]/sum(coun$V1)
  }
  coun$cell<-cell_type
  coun$Stat<-rownames(coun)
  data<-rbind(data,coun)
}
dim(all_data)
sum(data$V1)
dim(data)
# 查看结果

setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/细胞stat比例/补充")
write.table(data,"cell_stat_assintment训练.txt",row.names = FALSE)

###绘制堆叠箱线图
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
color<-c("S01" = "#D4256B", "S02" = "#3B95AE", "S03" = "#5CA84B", "S04" = "#754B9A", "S05" = "#FF8000")
p = ggplot( data, aes( x = cell, weight = percent, fill = Stat))+
  geom_bar( position = "stack")+scale_fill_manual( values = color)+labs(y = "Percent")+  theme_classic()+coord_flip()
ggsave("比例.pdf",plot = p,height=5,width=3.5)