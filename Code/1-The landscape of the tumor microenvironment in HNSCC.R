rm(list=ls())
library(Seurat)
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot)
library(ggpubr);library(beeswarm);library(presto);
library(GSEABase)
sc<-readRDS("/groups/g900008/home/zhangnihui/新课题/实验/12.13修改数据处理标准sample设置为HPV/训练集/recluster.rds")
scobj<-sc
#combined_data<-sc@meta.data

colors <- c(
  "Epithelial_cells" = "#699ECA",
  "Dendritic_cells" = "#FF8C00",
  "Fibroblasts" = "#F898CB",
  "Mono/Mac" = "#4DAF4A",
  "CD8_T_cells" = "#D65190",
  "Neutrophils" = "#731A73",
  "CD4_T_cells" = "#FFCB5B",
  "Endothelial_cells" = "#E87B1E",
  "B_cells" = "#0076B9",
  "NK_cells" = "#3D506A",
  "Mast_cells" = "#0099B2",
  "Plasma_cells" = "#F7B799"
)
col_sankey<-c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 
              'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono/Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF")
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/图片1")

##HPV分布
DimPlot(scobj,reduction = "umap",label=F,group.by="recluster")+scale_color_manual(values = col_sankey)
ggsave("recluster-umap.pdf",width=6,height=5)
head(scobj@meta.data)
colnames(scobj@meta.data)[4]<-"HPV_Status"
#saveRDS(scobj,"recluster2.rds")
colors1 <- c(
  "HPV+" = "#F8807B",
  "HPV-" = "#FFCB4E"
)
DimPlot(scobj,reduction = "umap",label=F,group.by="HPV_Status")+scale_color_manual(values = colors1)
ggsave("HPV-umap.pdf",width=5,height=5)



##特征基因表达
gene <- c('MS4A1','CD79A', # B_cells
          'MZB1','TNFRSF17','SDC1', # Plasma_cells
          'CD8A','CD8B', # CD8_T_cells
          'CD4', # CD4_T_cells
          'KLRF1', # NK_cells
          'AIF1','CD14', # Mono/Mac
          'S100A8','S100A12','CSF3R', # Neutrophils
          'CD1C','LY75','IRF4', # Dendritic_cells
          'CPA3','TPSB2','MS4A2', # Mast.cells
          'LUM','COL1A1','COL1A2', # Fibroblasts
          'VWF','ENG', # Endothelial_cells
          'KRT18','KRT5','KRT6A' # Epithelial_cells
          )
		
limiy<-c("B_cells","Plasma_cells","CD8_T_cells","CD4_T_cells","NK_cells","Mono/Mac","Neutrophils",
         "Dendritic_cells","Mast_cells","Fibroblasts","Endothelial_cells","Epithelial_cells") 
Idents(scobj)<-scobj$recluster
p<-DotPlot(scobj, features = gene)+
  scale_color_gradientn(colours = c('white','#E42320'))+#, limits = c(0, 2.5),
                        #na.value = "white")+
  #scale_x_discrete(limits = limix)+
  scale_y_discrete(limits = limiy)+ 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))+
  theme(axis.text.y = element_text(size = 12))+
  scale_size(range = c(1, 10))
ggsave("all_marker_qp.pdf",plot = p,height=5,width=12)


mytheme <- theme_bw() + 
  theme(plot.title=element_text(size=rel(2),hjust=0.5),
        axis.title=element_text(size=rel(1.2)),
        axis.text.x = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.5),color = 'black'),
        axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.5),color = 'black'),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        panel.border=element_rect(color="black",size=1),
        axis.line=element_line(color="black",size=0.5))


##细胞比例分布变化-堆叠柱状图
seucount<-scobj
N<-seucount@meta.data[which(seucount$HPV_Status=="HPV-"),]

P<-seucount@meta.data[which(seucount$HPV_Status=="HPV+"),]


#
coun<-as.data.frame(as.matrix(table(N$recluster),,2))
coun
for(i in 1:nrow(coun)){
  coun$percent[i]<-coun[i,1]/sum(coun$V1)
}
coun
#coun$ident<-"N01-nomal"
coun$sample<-"HPV-"
nomal<-coun

coun<-as.data.frame(as.matrix(table(P$recluster),,2))
coun
for(i in 1:nrow(coun)){
  coun$percent[i]<-coun[i,1]/sum(coun$V1)
}
coun
#coun$ident<-"N01-nomal"
coun$sample<-"HPV+"
patient<-coun

nomal$cell<-rownames(nomal)

patient$cell<-rownames(patient)

data<-cbind(nomal,patient)
data
data<-rbind(nomal,patient)
colors1 <- c(
  "HPV+" = "#F8807B",
  "HPV-" = "#FFCB4E"
)
data$sample<-factor(data$sample, level=c("HPV-","HPV+"))
p = ggplot( data, aes( x = sample, weight = percent, fill = cell))+
  geom_bar( position = "stack")+scale_fill_manual( values = col_sankey)+labs(y = "Percent")+mytheme
p
ggsave("比例.pdf",plot = p,height=6,width=4)


##细胞比例分布变化-箱线图
result<-c()
a<-unique(unique(N$orig.ident))
for(i in 1:length(a)){
  pri1<-N[which(N$orig.ident==a[i]),]
  number<-nrow(pri1)
  result1<-table(pri1$recluster)/number
  result1<-cbind(as.numeric(result1),names(result1),"HPV-")
  result<-rbind(result,result1)
}

a<-unique(unique(P$orig.ident))
for(i in 1:length(a)){
  pri1<-P[which(P$orig.ident==a[i]),]
  number<-nrow(pri1)
  result1<-table(pri1$recluster)/number
  result1<-cbind(as.numeric(result1),names(result1),"HPV+")
  result<-rbind(result,result1)
}
colnames(result)<-c("p","celltype","sample")
#colors<-c("#62B197","#E18E6D")

p <- ggplot(result, aes(x = celltype, y = as.numeric(p), fill = sample)) +
  geom_boxplot() +
  ylab("p") +
  xlab("cell type") +
  scale_fill_manual(values = c("Normal" = "blue", "Mild" = "green", "Severe" = "red")) +
  scale_fill_manual(values = colors1) +  # 应用之前定义的颜色向量
  theme(legend.title = element_blank()) +
  ylab("Proportion of cells") +  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # 设置X轴标签倾斜45度
library(ggpubr)

p <- p + stat_compare_means(label = "p.format",method = "wilcox.test")
ggsave("box-cellp-N-P.pdf",plot = p,height=4,width=12)



###平均值-细胞比例分布变化-箱线图
library(dplyr)
result_df <- as.data.frame(result)

# 分组并计算平均值
result_summary <- result_df %>%
  group_by(celltype, sample) %>%
  summarize(mean_p = mean(as.numeric(p)))

# 查看结果
print(result_summary)

p1<-ggplot(result_summary, aes(x = sample, y = mean_p, group = sample, color = celltype)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_line(aes(group = celltype), position = position_dodge(width = 0.7), size = 1) +
  scale_color_manual(values = col_sankey) +
  theme_classic() +
  labs(title = "Mean p-value by Sample and Cell Type", x = "Sample", y = "Mean p-value") +
  theme(legend.title = element_blank())
ggsave("box.pdf",plot = p1,height=3,width=3.5)
