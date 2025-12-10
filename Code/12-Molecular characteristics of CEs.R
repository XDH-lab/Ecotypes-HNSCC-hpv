#生存分析：
rm(list=ls())

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/bulk合并4/预后/ECO预后")

eco<-read.table("/lkn_lab/znh/12月课题实验/实验/验证集/bulk合并4/结果/exp_laster/Ecotypes/ecotype_assignment.txt")
head(eco)
lnc<-eco[,c(1,2,7,8)]

head(lnc)
table(lnc$MaxEcotype)
lnc<-na.omit(lnc)
colnames(lnc)[c(3,4)]<-c("fustat","futime")

is.numeric(lnc$futime)
is.numeric(lnc$fustat)

ls<-lnc
head(lnc)

library(survival)
survival_object <- Surv(lnc$futime, lnc$fustat)
fit <- survfit(survival_object ~ lnc$MaxEcotype)
library(ggplot2)

library(survminer)

setwd("/lkn_lab/znh/12月课题实验/实验/验证集/bulk合并4/预后/ECO预后/颜色")
col_sankey<-c('B_cells'="#E6AB02", 'CD4_T_cells'="#E7298A", 'CD8_T_cells'="#7570B3", 
              'Dendritic_cells'="#D95F02",'Endothelial_cells'="#1B9E77",'Epithelial_cells'="#CAB2D6",
              'Fibroblasts'="#FDBF6F",'Mast_cells'="#F7AEB3",
              'Mono_Mac'="#B2DF8A",'Neutrophils'="#A6CEE3",'NK_cells'="#999999",'Plasma_cells'="#F780BF",
              'E1'="#D6372E",'E2'="#FADD4B",'E3'="#70B460",'E4'="#E690C1",
              'E5'="#985EA8",'E6'="#A3A3A3")
ggsurvplot(fit, data = lnc, conf.int = FALSE,  # 关闭置信区间阴影
                     xlab = "Time", ylab = "Survival probability",
                     pval = TRUE,
					 palette = c("#D6372E", "#FADD4B","#70B460","#E690C1", "#985EA8","#A3A3A3"))
ggsave("surviva_eco.pdf")



#生态型功能变化
seucount<-readRDS("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/转录因子分析/eco_sc.rds")
library(ClusterGVis)
library(org.Hs.eg.db)
library(ComplexHeatmap)
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/ECO富集/修改")


###先找差异基因
data<-seucount
head(data@meta.data)
Idents(object = data) <- data$Ecotype
data_E1 <- subset(data, idents = "E1")
data_E2 <- subset(data, idents = "E2")
data_E3 <- subset(data, idents = "E3")
data_E4 <- subset(data, idents = "E4")
data_E5 <- subset(data, idents = "E5")
data_E6 <- subset(data, idents = "E6")
data_combined <- merge(data_E1, y = list(data_E2, data_E3, data_E4, data_E5,data_E6))
data<-data_combined
pbmc.markers.all <- Seurat::FindAllMarkers(data,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)
table(pbmc.markers.all$cluster)
# get top 50 genes
pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 50, wt = avg_log2FC)
# prepare data from seurat object

st.data <- prepareDataFromscRNA(object = data,
                                diffData = pbmc.markers,
                                showAverage = TRUE)

str(st.data)
ls<-str(st.data)

# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
#折线图
# add gene name
markGenes = pbmc.markers.all %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  ungroup()
markGenes<-markGenes$gene
# line plot
visCluster(object = st.data,
           plot.type = "line")

ggsave("line_plot.pdf")
# heatmap plot
pdf('sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:6))

dev.off()

Color1<-c("#D6372E","#FADD4B","#70B460","#E690C1",
           "#985EA8","#A3A3A3")
# add bar plot
pdf('marker_heatmap.pdf',height = 7,width = 12,onefile = F)

visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:6),
           go.col = rep(Color1,each = 5),
           add.bar = T)
dev.off()		   
		   
sample_order <- c("S05", "S04", "S03", "S02", "S01")
pdf('ls.pdf',height = 7,width = 12,onefile = F)

visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(5,1,2,3,4,6),
           go.col = rep(Color1,each = 5),
           add.bar = T)
dev.off()	