seucount<-readRDS("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/scdata/Epithelial_cells.rds")
library(ClusterGVis)
library(org.Hs.eg.db)
library(ComplexHeatmap)
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/scdata/富集/EP")

###先找差异基因
data<-seucount
Idents(object = data) <- data$stat
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
#data$stat<-factor(data$stat,levels=c("S01", "S02", "S03", "S04", "S05"))

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
markGenes = c(
  "SPRR3", "FCN1", "PRSS22", "TMPRSS11D", "PRSS27",
  "CES1", "ADH1C", "ALDH3A1", "AKR1C2", "AKR1B10",
  "IL1R2", "CXCL9", "PTN", "PRRX1", "CYP27A1",
  "NTS", "DAPL1", "FABP1", "DBI", "AQP3",
  "MAGEA4", "PCSK1N", "CXCL14", "CSAG1"
)
head(pbmc.markers)
markGenes<-intersect(pbmc.markers$gene, markGenes)
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
           cluster.order = c(1:5))

dev.off()

Color1 <- c("#8DD3C7", "#F1C900", "#938CC1", "#FB8072", "#629FC8")#五种细胞类型配色
Color1<-c("S01" = "#D4256B", "S02" = "#3B95AE", "S03" = "#5CA84B", "S04" = "#754B9A", "S05" = "#FF8000")
Color1<-c("#D4256B","#3B95AE","#5CA84B","#754B9A","#FF8000")

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
           cluster.order = c(1:5),
           go.col = rep(Color1,each = 5),
           add.bar = T)
dev.off()		   
		   
sample_order <- c("S05", "S04", "S03", "S02", "S01")


#######		   
seucount<-readRDS("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/scdata/NK_cells.rds")
library(ClusterGVis)
library(org.Hs.eg.db)
library(ComplexHeatmap)
setwd("/lkn_lab/znh/12月课题实验/实验/训练集/结果7进一步分析/scdata/富集/NK")


###先找差异基因
data<-seucount
s01_s02_data <- subset(data, subset = stat %in% c("S01", "S02","S03", "S04"))
head(s01_s02_data@meta.data)
Idents(object = data) <- data$stat
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

# check
str(st.data)

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
markGenes = c(
  "IL7R", "REL", "DUSP4", "ADGRE5", "CRTAM", 
  "TCF7", "CTSW", "KLRD1", "FGFBP2", "GZMH", 
  "FCGR3A", "S1PR5", "ITGB2", "CST7", "PRF1", 
  "NKG7", "CD8A", "GZMA", "GZMB", "RPS10", 
  "PSMA2", "NDUFB8", "RPS17", "TMSB4X", "CXCR6", 
  "HAVCR2", "KLRB1", "ITGA1", "PLPP1", "ITGAE", 
  "LDLRAD4", "LAYN", "PTPRCAP", "TIGIT"
)
head(pbmc.markers)
markGenes<-intersect(pbmc.markers$gene, markGenes)
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
           cluster.order = c(1:5))

dev.off()

Color1<-c("#D4256B","#3B95AE","#5CA84B","#754B9A")
# add bar plot
pdf('marker_heatmap.pdf',height = 6,width = 12,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:4),
           go.col = rep(Color1,each = 5),
           add.bar = T)
dev.off()